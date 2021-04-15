/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DDpackage_UNIQUETABLE_HPP
#define DDpackage_UNIQUETABLE_HPP

#include "ComplexNumbers.hpp"
#include "Definitions.hpp"
#include "Edge.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <iostream>
#include <limits>
#include <numeric>
#include <stack>
#include <stdexcept>
#include <vector>

namespace dd {

    /// Data structure for providing and uniquely storing DD nodes
    /// \tparam Node class of nodes to provide/store
    /// \tparam NBUCKET number of hash buckets to use (has to be a power of two)
    /// \tparam INITIAL_ALLOCATION_SIZE number if nodes initially allocated
    /// \tparam GROWTH_PERCENTAGE percentage that the allocations' size shall grow over time
    /// \tparam INITIAL_GC_LIMIT number of nodes initially used as garbage collection threshold
    /// \tparam GC_INCREMENT absolute number of nodes to increase the garbage collection threshold after garbage collection has been performed
    template<class Node, std::size_t NBUCKET = 32768, std::size_t INITIAL_ALLOCATION_SIZE = 2048, std::size_t GROWTH_FACTOR = 2, std::size_t INITIAL_GC_LIMIT = 250000, std::size_t GC_INCREMENT = 0>
    class UniqueTable {
    public:
        explicit UniqueTable(std::size_t nvars):
            nvars(nvars), chunkID(0), allocationSize(INITIAL_ALLOCATION_SIZE), gcLimit(INITIAL_GC_LIMIT) {
            // allocate first chunk of nodes
            chunks.emplace_back(allocationSize);
            allocations += allocationSize;
            allocationSize *= GROWTH_FACTOR;
            chunkIt    = chunks[0].begin();
            chunkEndIt = chunks[0].end();
        }

        ~UniqueTable() = default;

        static constexpr std::size_t MASK = NBUCKET - 1;

        void resize(std::size_t nq) {
            nvars = nq;
            tables.resize(nq);
            // TODO: if the new size is smaller than the old one we might have to release the unique table entries for the superfluous variables
            active.resize(nq);
            activeNodeCount = std::accumulate(active.begin(), active.end(), 0UL);
        }

        static std::size_t hash(const Node* p) {
            std::uintptr_t key = 0;
            for (std::size_t i = 0; i < p->e.size(); ++i) {
                key += ((reinterpret_cast<std::uintptr_t>(p->e[i].p) >> i) +
                        (reinterpret_cast<std::uintptr_t>(p->e[i].w.r) >> i) +
                        (reinterpret_cast<std::uintptr_t>(p->e[i].w.i) >> (i + 1))) &
                       MASK;
                key &= MASK;
            }
            return key;
        }

        // access functions
        [[nodiscard]] std::size_t getNodeCount() const { return nodeCount; }
        [[nodiscard]] std::size_t getPeakNodeCount() const { return peakNodeCount; }
        [[nodiscard]] std::size_t getAllocations() const { return allocations; }
        [[nodiscard]] float       getGrowthFactor() const { return GROWTH_FACTOR; }
        [[nodiscard]] const auto& getTables() const { return tables; }

        // lookup a node in the unique table for the appropriate variable; insert it, if it has not been found
        // NOTE: reference counting is to be adjusted by function invoking the table lookup and only normalized nodes shall be stored.
        Edge<Node> lookup(const Edge<Node>& e, bool keepNode = false) {
            // there are unique terminal nodes
            if (e.isTerminal())
                return e;

            lookups++;
            const auto key = hash(e.p);
            const auto v   = e.p->v;

            // successors of a node shall either have successive variable numbers or be terminals
            for ([[maybe_unused]] const auto& edge: e.p->e)
                assert(edge.p->v == v - 1 || edge.isTerminal());

            auto& bucket = tables[v][key];
            auto  it     = bucket.begin();
            while (it != bucket.end()) {
                if (e.p->e == (*it)->e) {
                    // Match found
                    if (e.p != (*it) && !keepNode) {
                        // put node pointed to by e.p on available chain
                        returnNode(e.p);
                    }
                    hits++;

                    // variables should stay the same
                    assert((*it)->v == e.p->v);

                    // successors of a node shall either have successive variable numbers or be terminals
                    for ([[maybe_unused]] const auto& edge: e.p->e)
                        assert(edge.p->v == v - 1 || edge.isTerminal());

                    return {(*it), e.w};
                }
                collisions++;
                ++it;
            }

            // node was not found -> add it to front of unique table bucket
            bucket.push_front(e.p);
            nodeCount++;
            peakNodeCount = std::max(peakNodeCount, nodeCount);

            return e;
        }

        [[nodiscard]] Node* getNode() {
            // a node is available on the stack
            if (!available.empty()) {
                auto p = available.top();
                available.pop();
                // returned nodes could have a ref count != 0
                p->ref = 0;
                return p;
            }

            // new chunk has to be allocated
            if (chunkIt == chunkEndIt) {
                chunks.emplace_back(allocationSize);
                allocations += allocationSize;
                allocationSize *= GROWTH_FACTOR;
                chunkID++;
                chunkIt    = chunks[chunkID].begin();
                chunkEndIt = chunks[chunkID].end();
            }

            auto p = &(*chunkIt);
            ++chunkIt;
            return p;
        }

        void returnNode(Node* p) {
            available.push(p);
        }

        // increment reference counter for node e points to
        // and recursively increment reference counter for
        // each child if this is the first reference
        void incRef(const Edge<Node>& e) {
            dd::ComplexNumbers::incRef(e.w);
            if (e.isTerminal())
                return;

            if (e.p->ref == std::numeric_limits<decltype(e.p->ref)>::max()) {
                std::clog << "[WARN] MAXREFCNT reached for p=" << reinterpret_cast<std::uintptr_t>(e.p) << ". Node will never be collected." << std::endl;
                return;
            }

            e.p->ref++;

            if (e.p->ref == 1) {
                for (const auto& edge: e.p->e) {
                    if (edge.p != nullptr) {
                        incRef(edge);
                    }
                }
                active[e.p->v]++;
                activeNodeCount++;
                maxActive = std::max(maxActive, activeNodeCount);
            }
        }

        // decrement reference counter for node e points to
        // and recursively decrement reference counter for
        // each child if this is the last reference
        void decRef(const Edge<Node>& e) {
            dd::ComplexNumbers::decRef(e.w);
            if (e.isTerminal()) return;
            if (e.p->ref == std::numeric_limits<decltype(e.p->ref)>::max()) return;

            if (e.p->ref == 0) {
                throw std::runtime_error("In decref: ref==0 before decref\n");
            }

            e.p->ref--;

            if (e.p->ref == 0) {
                for (const auto& edge: e.p->e) {
                    if (edge.p != nullptr) {
                        decRef(edge);
                    }
                }
                active[e.p->v]--;
                activeNodeCount--;
            }
        }

        std::size_t garbageCollect(bool force = false) {
            gcCalls++;
            if (!force && nodeCount < gcLimit)
                return 0;

            gcRuns++;
            std::size_t collected = 0;
            std::size_t remaining = 0;
            for (auto& table: tables) {
                for (auto& bucket: table) {
                    auto it     = bucket.begin();
                    auto lastit = bucket.before_begin();
                    while (it != bucket.end()) {
                        if ((*it)->ref == 0) {
                            assert(!Node::isTerminal(*it));

                            auto node = (*it);
                            bucket.erase_after(lastit); // erases the element at `it`
                            returnNode(node);
                            if (lastit == bucket.before_begin()) {
                                // first entry was removed
                                it = bucket.begin();
                            } else {
                                // entry in middle of list was removed
                                it = ++lastit;
                            }
                            collected++;
                        } else {
                            ++it;
                            ++lastit;
                            remaining++;
                        }
                    }
                }
            }
            gcLimit += GC_INCREMENT;
            nodeCount = remaining;
            return collected;
        }

        void clear() {
            // clear unique table buckets
            for (auto& table: tables) {
                for (auto& bucket: table) {
                    bucket.clear();
                }
            }
            // clear available stack
            while (!available.empty())
                available.pop();

            // release memory of all but the first chunk TODO: it could be desirable to keep the memory
            while (chunkID > 0) {
                chunks.pop_back();
                chunkID--;
            }
            // restore initial chunk setting
            chunkIt        = chunks[0].begin();
            chunkEndIt     = chunks[0].end();
            allocationSize = INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR;
            allocations    = INITIAL_ALLOCATION_SIZE;

            nodeCount     = 0;
            peakNodeCount = 0;

            collisions = 0;
            hits       = 0;
            lookups    = 0;

            std::fill(active.begin(), active.end(), 0);
            activeNodeCount = 0;
            maxActive       = 0;

            gcCalls = 0;
            gcRuns  = 0;
            gcLimit = INITIAL_GC_LIMIT;
        };

        void print() {
            Qubit q = nvars - 1;
            for (auto it = tables.rbegin(); it != tables.rend(); ++it) {
                auto& table = *it;
                std::cout << "\t" << static_cast<std::size_t>(q) << ":"
                          << "\n";
                for (std::size_t key = 0; key < table.size(); ++key) {
                    auto& bucket = table[key];
                    if (!bucket.empty())
                        std::cout << key << ": ";

                    for (const auto& node: bucket)
                        std::cout << "\t\t" << reinterpret_cast<std::uintptr_t>(node) << " " << node->ref << "\t";

                    if (!bucket.empty())
                        std::cout << "\n";
                }
                --q;
            }
        }

        void printActive() {
            std::cout << "#printActive: " << activeNodeCount << ", ";
            for (const auto& a: active)
                std::cout << a << " ";
            std::cout << "\n";
        }

        [[nodiscard]] fp hitRatio() const { return static_cast<fp>(hits) / lookups; }
        [[nodiscard]] fp colRatio() const { return static_cast<fp>(collisions) / lookups; }

        [[nodiscard]] std::size_t getActiveNodeCount() const {
            return activeNodeCount;
        }
        [[nodiscard]] std::size_t getActiveNodeCount(Qubit var) { return active.at(var); }

        std::ostream& printStatistics(std::ostream& os = std::cout) {
            os << "hits: " << hits << ", collisions: " << collisions << ", looks: " << lookups << ", hitRatio: " << hitRatio() << ", colRatio: " << colRatio() << ", gc calls: " << gcCalls << ", gc runs: " << gcRuns << "\n";
            return os;
        }

    private:
        using NodeBucket = std::forward_list<Node*>;
        using Table      = std::array<NodeBucket, NBUCKET>;

        // unique tables (one per input variable)
        std::size_t        nvars = 0;
        std::vector<Table> tables{nvars};

        std::stack<Node*>                    available{};
        std::vector<std::vector<Node>>       chunks{};
        std::size_t                          chunkID;
        typename std::vector<Node>::iterator chunkIt;
        typename std::vector<Node>::iterator chunkEndIt;
        std::size_t                          allocationSize;

        std::size_t allocations   = 0;
        std::size_t nodeCount     = 0;
        std::size_t peakNodeCount = 0;

        // unique table lookup statistics
        std::size_t collisions = 0;
        std::size_t hits       = 0;
        std::size_t lookups    = 0;

        // (max) active nodes
        // number of active vector nodes for each variable
        std::vector<std::size_t> active{nvars};
        std::size_t              activeNodeCount = 0;
        std::size_t              maxActive       = 0;

        // garbage collection
        std::size_t gcCalls = 0;
        std::size_t gcRuns  = 0;
        std::size_t gcLimit = 250000;
    };

} // namespace dd

#endif //DDpackage_UNIQUETABLE_HPP