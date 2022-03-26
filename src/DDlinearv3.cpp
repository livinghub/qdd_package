#include "DDpackage.h"


namespace dd {
    Edge Package::linearSifting(Edge in, std::map<unsigned short, unsigned short>& varMap) {
        const auto n = static_cast<short>(in.p->v + 1); //变量个数

		std::vector<bool> free(n, true); //记录数组
		std::map<unsigned short, unsigned short> invVarMap{}; //DD qubit（变量） 到 电路 qubit（变量） 的映射
		for (const auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边

        unsigned int total_min = size(in);
        unsigned int index_min = 0, index_mid=0; //开区间索引

        short pos = -1; //circuit variable (index)
        std::vector<short> moves;
        for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            // unsigned long min = size(in);
            unsigned long max = 0;

            // std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的变量为处理状态
            short optimalPos = pos; 
            short originalPos = pos;

            if(pos == n-1) {
                // linear sifting to bottom
                while (pos > 0) {
                    exchangeBaseCase(pos, in, varMap);
                    auto ex_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    linearInPlace(pos, in, varMap); //再pos，pos-1执行L.T.
                    auto lt_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    if(ex_size <= lt_size){
                        linearInPlace(pos, in, varMap);
                        assert(is_locally_consistent_dd(in));
                        moves.push_back(pos);

                        if(ex_size <= total_min)  {
                            total_min = ex_size;
                            index_min = moves.size();
                        }
                    } else {
                        moves.push_back(pos);
                        moves.push_back(-pos);

                        if(lt_size <= total_min)  {
                            total_min = lt_size;
                            index_min = moves.size();
                        }
                    }

                    //std::clog << "↓" << in_size << " ";
                    --pos;
                    
                } //到这里被选中的变量走到了最下面
            } else if(pos == 0) {
                //linear sifting top
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto ex_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    linearInPlace(pos+1, in, varMap); //再pos，pos-1执行L.T.
                    auto lt_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    if(ex_size <= lt_size){
                        linearInPlace(pos+1, in, varMap);
                        assert(is_locally_consistent_dd(in));
                        moves.push_back(pos+1);

                        if(ex_size <= total_min)  {
                            total_min = ex_size;
                            index_min = moves.size();
                        }

                    } else {
                        moves.push_back((pos+1));
                        moves.push_back(-(pos+1));

                        if(lt_size <= total_min)  {
                            total_min = lt_size;
                            index_min = moves.size();
                        }
                    }
                    //std::clog << "↑" << in_size << " ";
                    ++pos;
                }
            } else if (pos < n / 2) {  // variable is in lower half -> sifting to bottom first
                // linear sifting to bottom
                while (pos > 0) {
                    exchangeBaseCase(pos, in, varMap);
                    auto ex_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    linearInPlace(pos, in, varMap); //再pos，pos-1执行L.T.
                    auto lt_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    if(ex_size <= lt_size){
                        linearInPlace(pos, in, varMap);
                        assert(is_locally_consistent_dd(in));
                        moves.push_back(pos);

                        if(ex_size <= total_min)  {
                            total_min = ex_size;
                            index_min = moves.size();
                        }
                    } else {
                        moves.push_back(pos);
                        moves.push_back(-pos);

                        if(lt_size <= total_min)  {
                            total_min = lt_size;
                            index_min = moves.size();
                        }
                    }

                    //std::clog << "↓" << in_size << " ";
                    --pos;
                    
                } //到这里被选中的变量走到了最下面

                //还原到初始状态
                index_mid = moves.size();
                if(index_mid > 0) {
                    for(unsigned int i=index_mid-1; i>=0; --i) {
                        if(moves[i] > 0) {
                            exchangeBaseCase(moves[i], in, varMap);
                            assert(is_locally_consistent_dd(in));
                        } else {
                            linearInPlace(-moves[i], in, varMap);
                            assert(is_locally_consistent_dd(in));
                        }
                    }
                }

                //linear sifting top
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto ex_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    linearInPlace(pos+1, in, varMap); //再pos，pos-1执行L.T.
                    auto lt_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    if(ex_size <= lt_size){
                        linearInPlace(pos+1, in, varMap);
                        assert(is_locally_consistent_dd(in));
                        moves.push_back(pos+1);

                        if(ex_size <= total_min)  {
                            total_min = ex_size;
                            index_min = moves.size();
                        }

                    } else {
                        moves.push_back((pos+1));
                        moves.push_back(-(pos+1));

                        if(lt_size <= total_min)  {
                            total_min = lt_size;
                            index_min = moves.size();
                        }
                    }
                    //std::clog << "↑" << in_size << " ";
                    ++pos;
                }

                
            } else {
                 //linear sifting top
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto ex_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    linearInPlace(pos+1, in, varMap); //再pos，pos-1执行L.T.
                    auto lt_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    if(ex_size <= lt_size){
                        linearInPlace(pos+1, in, varMap);
                        assert(is_locally_consistent_dd(in));
                        moves.push_back(pos+1);

                        if(ex_size <= total_min)  {
                            total_min = ex_size;
                            index_min = moves.size();
                        }

                    } else {
                        moves.push_back((pos+1));
                        moves.push_back(-(pos+1));

                        if(lt_size <= total_min)  {
                            total_min = lt_size;
                            index_min = moves.size();
                        }
                    }
                    //std::clog << "↑" << in_size << " ";
                    ++pos;
                }
                //还原到初始状态
                index_mid = moves.size();
                if(index_mid > 0) {
                    for(unsigned int i=index_mid-1; i>=0; --i) {
                        if(moves[i] > 0) {
                            exchangeBaseCase(moves[i], in, varMap);
                            assert(is_locally_consistent_dd(in));
                        } else {
                            linearInPlace(-moves[i], in, varMap);
                            assert(is_locally_consistent_dd(in));
                        }
                    }
                }
                // linear sifting to bottom
                while (pos > 0) {
                    exchangeBaseCase(pos, in, varMap);
                    auto ex_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    linearInPlace(pos, in, varMap); //再pos，pos-1执行L.T.
                    auto lt_size = size(in);
                    assert(is_locally_consistent_dd(in));

                    if(ex_size <= lt_size){
                        linearInPlace(pos, in, varMap);
                        assert(is_locally_consistent_dd(in));
                        moves.push_back(pos);

                        if(ex_size <= total_min)  {
                            total_min = ex_size;
                            index_min = moves.size();
                        }
                    } else {
                        moves.push_back(pos);
                        moves.push_back(-pos);

                        if(lt_size <= total_min)  {
                            total_min = lt_size;
                            index_min = moves.size();
                        }
                    }

                    //std::clog << "↓" << in_size << " ";
                    --pos;
                    
                } //到这里被选中的变量走到了最下面
            }

            // linear sifting to optimal position
            if(index_min >= index_mid) {
                for(int i=moves.size()-1; i>=index_min; --i) {
                    if(moves[i] > 0) {
                        exchangeBaseCase(moves[i], in, varMap);
                        assert(is_locally_consistent_dd(in));
                    } else {
                        linearInPlace(-moves[i], in, varMap);
                        assert(is_locally_consistent_dd(in));
                    }
                }
            } else {
                for(int i=moves.size()-1; i>=index_mid; --i) {
                    if(moves[i] > 0) {
                        exchangeBaseCase(moves[i], in, varMap);
                        assert(is_locally_consistent_dd(in));
                    } else {
                        linearInPlace(-moves[i], in, varMap);
                        assert(is_locally_consistent_dd(in));
                    }
                }
                for(int i=0; i<index_min; ++i) {
                    if(moves[i] > 0) {
                        exchangeBaseCase(moves[i], in, varMap);
                        assert(is_locally_consistent_dd(in));
                    } else {
                        linearInPlace(-moves[i], in, varMap);
                        assert(is_locally_consistent_dd(in));
                    }
                }
            }

            initComputeTable();

                        // there are nodes which need to renormalized
            if (unnormalizedNodes > 0) {
                std::clog << "{" << unnormalizedNodes << "} ";
                auto oldroot = root;
                root = renormalize(root);
                decRef(oldroot);
                incRef(root);
                in.p = root.p;
                in.w = root.w;
                if (unnormalizedNodes > 0) {
                    throw std::runtime_error("Renormalization failed. " + std::to_string(unnormalizedNodes) + " unnormalized nodes remaining.");
                }
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root); //标记
            recomputeMatrixProperties(root);
        }
        return in;
    }

}