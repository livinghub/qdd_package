/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_COMPLEXCACHE_HPP
#define DD_PACKAGE_COMPLEXCACHE_HPP

#include "Complex.hpp"
#include "ComplexTable.hpp"

#include <cassert>
#include <cstddef>
#include <vector>

namespace dd {

    template<std::size_t INITIAL_ALLOCATION_SIZE = 2048, std::size_t GROWTH_FACTOR = 2>
    class ComplexCache {
        using Entry = ComplexTable<>::Entry;

    public:
        //构造函数,初始化
        ComplexCache():
            chunkID(0), allocationSize(INITIAL_ALLOCATION_SIZE) {
            // allocate first chunk of cache entries
            chunks.emplace_back(allocationSize);
            allocations += allocationSize;
            allocationSize *= GROWTH_FACTOR;
            chunkIt    = chunks[0].begin();
            chunkEndIt = chunks[0].end();
        }

        ~ComplexCache() = default;

        // access functions
        [[nodiscard]] std::size_t getCount() const { return count; }
        [[nodiscard]] std::size_t getPeakCount() const { return peakCount; }
        [[nodiscard]] std::size_t getAllocations() const { return allocations; }
        [[nodiscard]] std::size_t getGrowthFactor() const { return GROWTH_FACTOR; }

        //
        [[nodiscard]] Complex getCachedComplex() {
            // an entry is available on the stack
            if (available != nullptr) {
                assert(available->next != nullptr);
                auto entry = Complex{available, available->next}; //定义一个新记录
                available  = entry.i->next;
                count += 2;
                return entry;
            }

            // new chunk has to be allocated
            if (chunkIt == chunkEndIt) { //内存区用完了
                chunks.emplace_back(allocationSize); //在尾部新增一个序号为allocationSize的内存区
                allocations += allocationSize; //更新总的内存分配空间
                allocationSize *= GROWTH_FACTOR; //为下一个内存区作准备
                chunkID++; //增加内存区计数
                chunkIt    = chunks[chunkID].begin(); //内存区内第一个元素(内存块)
                chunkEndIt = chunks[chunkID].end(); //内存区内最后一个元素(内存块)
            }

            Complex c{};
            c.r = &(*chunkIt); //实部虚部各占一个内存块
            ++chunkIt;
            c.i = &(*chunkIt);
            ++chunkIt;
            count += 2;
            return c;
        }

        [[nodiscard]] Complex getTemporaryComplex() {
            // an entry is available on the stack
            if (available != nullptr) {
                assert(available->next != nullptr);
                return {available, available->next};
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
            return {&(*chunkIt), &(*(chunkIt + 1))};
        }

        //退出一个cache
        void returnToCache(Complex& c) {
            assert(count >= 2);
            assert(c != Complex::zero);
            assert(c != Complex::one);
            assert(c.r->refCount == 0);
            assert(c.i->refCount == 0);
            c.i->next = available;
            c.r->next = c.i;
            available = c.r;
            count -= 2;
        }

        void clear() {
            // clear available stack
            available = nullptr;

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

            count     = 0;
            peakCount = 0;
        };

    private:
        Entry*                                available{};
        std::vector<std::vector<Entry>>       chunks{};
        std::size_t                           chunkID;
        typename std::vector<Entry>::iterator chunkIt;
        typename std::vector<Entry>::iterator chunkEndIt;
        std::size_t                           allocationSize;

        std::size_t allocations = 0;
        std::size_t count       = 0;
        std::size_t peakCount   = 0;
    };
} // namespace dd

#endif //DD_PACKAGE_COMPLEXCACHE_HPP
