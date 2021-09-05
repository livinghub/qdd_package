#include "DDpackage.h"

namespace dd {

	void Package::printLTMap(std::map<unsigned short, unsigned short>& varMap)
	{
		unsigned short nqubits = varMap.size();
		for(int i=0; i<nqubits; ++i)
		{
			std::cout << i << "	:	";
			for(int j=0; j<nqubits; ++j)
			{
				if(xorMat[i][j] > 0)
					std::cout << varMap.at(j) << " ";
			}
			std::cout << std::endl;
		}
	}

	void Package::printLTMat(std::map<unsigned short, unsigned short>& varMap)
	{
		unsigned short nqubits = varMap.size();
		for(int i=0; i<nqubits; ++i)
		{
			for(int j=0; j<nqubits; ++j)
			{
				std::cout << Package::xorMat[varMap.at(i)][varMap.at(j)] << " ";
			}
			std::cout << std::endl;
		}
	}

    void Package::xorInit(std::map<unsigned short, unsigned short>& varMap) { //初始化lt矩阵
		const unsigned short varNum = varMap.size();
        for(int i=0; i<varNum; ++i) {
            Package::xorMat[i][i] = 1;
        }
    }

    void Package::xorLinear(unsigned short index) { //更新lt矩阵
        for(int i=0; i<MAXN; ++i)
            Package::xorMat[index][i] ^= Package::xorMat[index-1][i];
    }

	//更新lt矩阵，只维护一个lt矩阵
	void Package::xorLinear(unsigned short index, std::map<unsigned short, unsigned short>& varMap) { //更新lt矩阵
        const unsigned short varNum = varMap.size();
		std::pair<short, short> xorEntry(varMap[index], varMap[index-1]);
		//debug
		// for(const auto& Q: varMap) {
		// 	std::clog <<"\t" << Q.first << ": " << Q.second << std::endl;
		// }

		for(int i=0; i<varNum; ++i)
            xorMat[varMap[index]][i] ^= xorMat[varMap[index-1]][i];
			//xorMat[index][i] ^= xorMat[index-1][i];

		//记录xor的组合和顺序
		if(!LTpath.empty() && LTpath.back().first==xorEntry.first && LTpath.back().second==xorEntry.second) {
			LTpath.pop_back();
		} else {
			LTpath.emplace_back(xorEntry);
		}
    }


	void Package::linearInPlace(unsigned short i, Edge in, std::map<unsigned short, unsigned short>& varMap, bool updateLTMat) { //传入变量的索引和dd的根指针,和变量映射
		std::clog << "lt:" << i << ", ";

		// copy unique table from higher variable and empty it
		std::array<NodePtr, NBUCKET> table{};
		for (unsigned short bucket=0; bucket < NBUCKET; ++bucket) {
			table.at(bucket) = Unique[i][bucket]; //把i变量的唯一表拿出来
			Unique[i][bucket] = nullptr; //把属于i变量的唯一表清空
		} //到此，已经把i变量的唯一表拿出来了，并且清空了它

		initComputeTable(); //初始化计算表

		// iterate over all obtained nodes
		for (unsigned short bucket=0; bucket < NBUCKET; ++bucket) { //对i变量的每个哈希桶进行遍历
			NodePtr p = table[bucket]; //取出当前桶里的冲突链
			while (p != nullptr) { //遍历冲突链的每个结点
				NodePtr pnext = p->next;
				assert(p->v == i); //确认当前指针指向的结点就是属于i变量的
                assert(i-1 == p->e[0].p->v || isTerminal(p->e[0])); //确认这个结点的四个孩子结点都是属于i-1变量，或是它的孩子是终端结点
                assert(i-1 == p->e[1].p->v || isTerminal(p->e[1]));
                assert(i-1 == p->e[2].p->v || isTerminal(p->e[2]));
                assert(i-1 == p->e[3].p->v || isTerminal(p->e[3]));
				if (p->ref != 0) { //确认这个结点引用计数不是0
                    linearInPlace2(p, i, in); //传入该结点的指针，变量索引，DD指针
				}
                assert(p->v == i);
                assert(i-1 == p->e[0].p->v || isTerminal(p->e[0]));
                assert(i-1 == p->e[1].p->v || isTerminal(p->e[1]));
                assert(i-1 == p->e[2].p->v || isTerminal(p->e[2]));
                assert(i-1 == p->e[3].p->v || isTerminal(p->e[3]));
				p = pnext; //轮到冲突链的下一个结点
			}
		}

		linear_in_place++;
		if(updateLTMat) {
			xorLinear(i, varMap); //更新线性变换矩阵
		}
		
	}

	void Package::linearInPlace(unsigned short i, Edge in) { //传入变量的索引和dd的根指针
		std::clog << "lt:" << i << ", ";

		// linear_in_place++;
		// xorLinear(i); //更新线性变换矩阵

		// copy unique table from higher variable and empty it
		std::array<NodePtr, NBUCKET> table{};
		for (unsigned short bucket=0; bucket < NBUCKET; ++bucket) {
			table.at(bucket) = Unique[i][bucket]; //把i变量的唯一表拿出来
			Unique[i][bucket] = nullptr; //把属于i变量的唯一表清空
		} //到此，已经把i变量的唯一表拿出来了，并且清空了它

		initComputeTable(); //初始化计算表

		// iterate over all obtained nodes
		for (unsigned short bucket=0; bucket < NBUCKET; ++bucket) { //对i变量的每个哈希桶进行遍历
			NodePtr p = table[bucket]; //取出当前桶里的冲突链
			while (p != nullptr) { //遍历冲突链的每个结点
				NodePtr pnext = p->next;
				assert(p->v == i); //确认当前指针指向的结点就是属于i变量的
                assert(i-1 == p->e[0].p->v || isTerminal(p->e[0])); //确认这个结点的四个孩子结点都是属于i-1变量，或是它的孩子是终端结点
                assert(i-1 == p->e[1].p->v || isTerminal(p->e[1]));
                assert(i-1 == p->e[2].p->v || isTerminal(p->e[2]));
                assert(i-1 == p->e[3].p->v || isTerminal(p->e[3]));
				if (p->ref != 0) { //确认这个结点引用计数不是0
                    linearInPlace2(p, i, in); //传入该结点的指针，变量索引，DD指针
				}
                assert(p->v == i);
                assert(i-1 == p->e[0].p->v || isTerminal(p->e[0]));
                assert(i-1 == p->e[1].p->v || isTerminal(p->e[1]));
                assert(i-1 == p->e[2].p->v || isTerminal(p->e[2]));
                assert(i-1 == p->e[3].p->v || isTerminal(p->e[3]));
				p = pnext; //轮到冲突链的下一个结点
			}
		}
	}

	void Package::linearInPlace2(NodePtr p, unsigned short index, Edge in) {
		Edge t[NEDGE][NEDGE]{ }; //创建一个exchange用的矩阵
        //Edge LTt[NEDGE][NEDGE]{ }; //lt用到的矩阵
		assert(index > 0);
		assert(index == p->v); //确定传入结点是属于index变量
		//xorLinear(index);

		//for lt
		int x,y;
        // creating matrix T
		for (int i = 0; i < NEDGE; i++) { //遍历要处理结点的每条出边
			x = i;
			y = 0;
			for (int j = 0; j < NEDGE; j++) { //遍历它的孩子结点的每一条出边
				if (p->e[i].p->v == index - 1) { //如果要处理结点的孩子结点是属于它的下一个变量的（如本结点是i变量的，它的孩子结点是i-1变量的）
				    assert(!isTerminal(p->e[i]));

					t[x][y] = p->e[i].p->e[j]; //结点第i条边的孩子的第j条边放入矩阵中
					auto c = cn.getTempCachedComplex(); //拿到一个临时复数c
					CN::mul(c, p->e[i].p->e[j].w, p->e[i].w); //复数乘法，结点出边权值*它的孩子的出边权值=c
					if (p->e[i].p->normalizationFactor != CN::ONE) { //如果改结点的规范因子不是1，就把c乘上规范因子
						CN::mul(c, c, p->e[i].p->normalizationFactor);
					}
					t[x][y].w = cn.lookup(c); //给记录的边附上权值c
				} else if (isTerminal(p->e[i])) { //要处理结点的孩子结点是终端结点
					// edge pointing to a terminal
					t[x][y] = p->e[i]; //矩阵记录的就是结点的出边
                    assert(t[x][y].p->normalizationFactor == CN::ONE); //这样的话确定这条边的结点规范因子是1
				} else { //上面两种情况都不是的话，就出问题了，出错处理
				    debugnode(p);
				    std::stringstream hex_addr;
				    hex_addr << "0x" << std::hex << reinterpret_cast<std::uintptr_t>(p);
				    throw std::runtime_error("Edge " + std::to_string(i)
				        + " of " + hex_addr.str()
				        + " pointing to a skipped variable: "
				        + std::to_string(index) + " --> " + std::to_string(p->e[i].p->v));
				}
				//for lt
				if(i%2==0) {
					x = (x+1+NEDGE)%NEDGE;
					//y = (y+1+NEDGE)%NEDGE;
					++y;
				} else {
					x = (x-1+NEDGE)%NEDGE;
					//y = (y+1+NEDGE)%NEDGE;
					++y;
				}
			}
			
			
		} //到这里，要处理的结点和它的孩子结点全都处理完，并把结点的孩子结点的出边和联合权值放入矩阵
        //std::memcpy(t, LTt, sizeof(t)); //复制t矩阵给tmpt

        // for(int x=0, i, j; x<NEDGE; ++x) { //重排矩阵方便处理
        //     i=x, j=0;
        //     if(x%2 == 0)
        //         for(int y=0; y<NEDGE; ++y) LTt[y][x] = t[(NEDGE+i++)%NEDGE][j++];
        //     else 
        //         for(int y=0; y<NEDGE; ++y) LTt[y][x] = t[(NEDGE+i--)%NEDGE][j++];
        // }

        assert(is_locally_consistent_dd({p, CN::ZERO}));

        // creating new nodes and appending corresponding edges
		Edge newEdges[NEDGE]{ }; //创建要处理结点（P）的一组（4条）出边
		for (int x = 0; x < NEDGE; ++x) { //搞定p结点的四个孩子
			newEdges[x] = makeNonterminal(static_cast<short>(index - 1), t[x]); //规范新的index-1结点的四个孩子结点，输入index-1和index-1（4个）结点里面的第x条边
            incRef(newEdges[x]); //给新结点加引用计数
			assert(is_locally_consistent_dd(newEdges[x]));
		}//到此，变量索引为index-1的四个结点完成规范化

        for (dd::Edge& x : p->e) //把结点原来只想它的孩子的边撤掉
            decRef(x);
		// reuse p to build new top node
		assert(p->ref > 0);
        reuseNonterminal(static_cast<short>(index), newEdges, p, in); //对p（要处理的结点）重新规范化,传入index，刚刚生成新”第二层“结点，第一层结点指针，DD指针
		// p might be discarded at this point if nodes were substituted
	}



    ///linear sifting 
	std::tuple<Edge, unsigned int, unsigned int> Package::linearSifting(Edge in, std::map<unsigned short, unsigned short>& varMap) {
		const auto n = static_cast<short>(in.p->v + 1); //变量个数

		std::vector<bool> free(n, true); //记录数组
		std::map<unsigned short, unsigned short> invVarMap{}; //DD qubit（变量） 到 电路 qubit（变量） 的映射
		for (const auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边

		//std::clog << "  Start Sifting. n=" << std::setw(2) << n << " -- ";
//		for (auto &entry: varMap) {
		    //std::clog << entry.second << " ";
//		}
		//std::clog << "\n";

		unsigned int total_max = size(in); //DD的大小
		unsigned int total_min = total_max;

        short pos = -1;

		/**/
		
		xorInit(varMap);

		/**/


		for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;
			std::vector<bool> ltFlag(n, false);

            std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > (unsigned short)max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的变量为处理状态
            short optimalPos = pos; 
            short originalPos = pos;

			if (pos < n / 2) {  // variable is in lower half -> sifting to bottom first
                // sifting to bottom 向下尝试
                while (pos > 0) {
					std::clog << "向下尝试";
                    exchangeBaseCase(pos, in); //先对pos，pos-1执行swap
                    auto ex_size = size(in);
					//linearInPlace(pos, in); //再pos，pos-1执行L.T.
					linearInPlace(pos, in, varMap); //再pos，pos-1执行L.T.
					auto lt_size = size(in);

                    total_min = std::min({total_min, ex_size, lt_size}); //记录ls过程中产生的最大和最小size
                    total_max = std::max({total_max, ex_size, lt_size});

                    //std::clog << "↓" << ex_size << " ";
                    assert(is_locally_consistent_dd(in));
                    --pos; //变量位置下移一位

					if(ex_size <= lt_size){ //swap效果更好
						//** 抵消lt
						std::clog << "抵消lt";
						//linearInPlace(pos+1, in);
						linearInPlace(pos+1, in, varMap);
						if (ex_size < min) {
							min = ex_size;
							optimalPos = pos;
						}
					} else { //lt效果更好
						//** 记录这个步骤用了lt
						ltFlag.at(pos+1) = true;
						++valid_LT_Num;
						if (lt_size < min) {
							min = lt_size;
							optimalPos = pos;
						}
					}
                    
                } //到这里被选中的变量走到了最下面

				// 还原到初始状态 向上还原
				while (pos < originalPos)
				{
					std::clog << "向上还原";
					if(ltFlag.at(pos+1)) //如果用了lt
					{
						std::clog << "and LT";
						//linearInPlace(pos+1, in);
						linearInPlace(pos+1, in, varMap);
						--valid_LT_Num;
						ltFlag.at(pos+1) = false;
					}
					exchangeBaseCase(pos+1, in);
					//exchangeBaseCase(pos+1, in);
					auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
					++pos;
				}
				//std::clog << "+++++++++++" << std::endl;

                // sifting to top 向上尝试
                while (pos < n - 1) {
					std::clog << "向上尝试";
                    exchangeBaseCase(pos+1, in); //先执行swap
                    auto ex_size = size(in);
					//linearInPlace(pos+1, in);
					linearInPlace(pos+1, in, varMap);
					auto lt_size = size(in);
							

                    total_min = std::min({total_min, ex_size, lt_size}); //记录ls过程中产生的最大和最小size
                    total_max = std::max({total_max, ex_size, lt_size});
                    //std::clog << "↑" << ex_size << " ";
                    assert(is_locally_consistent_dd(in));
                    ++pos; //变量位置上移一位
					
					if(ex_size <= lt_size){ //swap效果更好
						//** 抵消lt
						std::clog << "抵消lt";
						//linearInPlace(pos, in);
						linearInPlace(pos, in, varMap);
						if (ex_size < min) {
							min = ex_size;
							optimalPos = pos;
						}
					} else { //lt效果更好
						//** 记录这个步骤用了lt
						ltFlag.at(pos) = true;
						++valid_LT_Num;
						if (lt_size < min) {
							min = lt_size;
							optimalPos = pos;
						}
					}
                } //到顶了

                //std::clog << "[" << min << "] ";

                // sifting to optimal position 向下还原 找到最佳位置
                while (pos > optimalPos) {
					std::clog << "向下还原";
					if(ltFlag.at(pos)) //如果到这个位置用了lt，那么需要把它抵消掉
					{ 
						std::clog << "and LT";
						//linearInPlace(pos, in);
						linearInPlace(pos, in, varMap);
						--valid_LT_Num;
						ltFlag.at(pos) = false;
					}
                    exchangeBaseCase(pos, in); //再次执行swap抵消操作
					//exchangeBaseCase(pos, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    //std::clog << "↓" << size(in) << " ";
                    assert(is_locally_consistent_dd(in));
                    --pos;
                }
            } else {  // variable is in upper half -> sifting to top first

                // sifting to top 向上尝试
                while (pos < n - 1) {
					std::clog << "向上尝试";
                    exchangeBaseCase(pos+1, in); //先执行swap
                    auto ex_size = size(in);
					//linearInPlace(pos+1, in); //再执行L.T.
					linearInPlace(pos+1, in, varMap);
					auto lt_size = size(in);

                    total_min = std::min({total_min, ex_size, lt_size}); //记录ls过程中产生的最大和最小size
                    total_max = std::max({total_max, ex_size, lt_size});
                    //std::clog << "↑" << ex_size << " ";
                    assert(is_locally_consistent_dd(in));
                    ++pos; //变量位置上移一位
					
					if(ex_size <= lt_size){ //swap效果更好
						//** 抵消lt
						std::clog << "抵消lt";
						//linearInPlace(pos, in);
						linearInPlace(pos, in, varMap);
						if (ex_size < min) {
							min = ex_size;
							optimalPos = pos;
						}
					} else { //lt效果更好
						//** 记录这个步骤用了lt
						ltFlag.at(pos) = true;
						++valid_LT_Num;
						if (lt_size < min) {
							min = lt_size;
							optimalPos = pos;
						}
					}
                }

				// 还原到初始状态 向下还原
				while (pos > originalPos)
				{
					std::clog << "向下还原";
					if(ltFlag.at(pos))
					{
						std::clog << "and LT";
						//linearInPlace(pos, in);
						linearInPlace(pos, in, varMap);
						--valid_LT_Num;
						ltFlag.at(pos) = false;
					}
					exchangeBaseCase(pos, in);
					//exchangeBaseCase(pos, in);
					auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
					--pos;
				}
				//std::clog << "-------------------" << std::endl;
                // sifting to bottom 向下尝试
                while (pos > 0) {
					std::clog << "向下尝试";
                    exchangeBaseCase(pos, in); //先执行swap
                    auto ex_size = size(in);
					//linearInPlace(pos, in); //再执行L.T.
					linearInPlace(pos, in, varMap);
					auto lt_size = size(in);

                    total_min = std::min({total_min, ex_size, lt_size}); //记录ls过程中产生的最大和最小size
                    total_max = std::max({total_max, ex_size, lt_size});

                    //std::clog << "↓" << ex_size << " ";
                    assert(is_locally_consistent_dd(in));
                    --pos; //变量位置下移一位

					if(ex_size <= lt_size){ //swap效果更好
						//** 抵消lt
						std::clog << "抵消lt";
						//linearInPlace(pos+1, in); 
						linearInPlace(pos+1, in, varMap);
						if (ex_size < min) {
							min = ex_size;
							optimalPos = pos;
						}
					} else { //lt效果更好
						//** 记录这个步骤用了lt
						ltFlag.at(pos+1) = true;
						++valid_LT_Num;
						if (lt_size < min) {
							min = lt_size;
							optimalPos = pos;
						}
					}
                    
                } //到这里被选中的变量走到了最下面

                //std::clog << "[" << min << "] ";

                // sifting to optimal position  向上还原
                while (pos < optimalPos) {
					std::clog << "向上还原";
					if(ltFlag.at(pos+1)) 
					{
						std::clog << "and LT";
						//linearInPlace(pos+1, in);
						linearInPlace(pos+1, in, varMap);
						--valid_LT_Num;
						ltFlag.at(pos+1) = false;
					}
                    exchangeBaseCase(pos+1, in);
					//exchangeBaseCase(pos+1, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    //std::clog << "↑" << size(in) << " ";
                    assert(is_locally_consistent_dd(in));
                    ++pos;
                }
            }

			initComputeTable();

			// there are nodes which need to renormalized
            if (unnormalizedNodes > 0) {
                //std::clog << "{" << unnormalizedNodes << "} ";
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

			// Adjusting varMap if position changed
            if (optimalPos != originalPos) {
                //std::clog << "| " << originalPos << "-->" << optimalPos << " (min=" << min << "; real size=" << size(in) << ")\n";
            } else {
                //std::clog << "| ##### (min=" << min << "; real size=" << size(in) << ")\n";
            }
			//** 这里可能还要修改 --- 改变varMap
			if (optimalPos > originalPos) { //向上到最佳位置
                auto tempVar = invVarMap[originalPos]; //暂存最佳位置对应的电路变量
                for (int j = originalPos; j < optimalPos; ++j) {
                    invVarMap[j] = invVarMap[j + 1]; //调整 电路变量（qubit）
                    varMap[invVarMap[j]] = j; //更新DD变量（qubit）
                }
                invVarMap[optimalPos] = tempVar; //在合适位置放入最佳位置对应的电路变量
                varMap[invVarMap[optimalPos]] = optimalPos; //在合适位置放入最佳位置对应的DD变量
            } else if (optimalPos < originalPos) {
                auto tempVar = invVarMap[originalPos];
                for (int j = originalPos; j > optimalPos; --j) {
                    invVarMap[j] = invVarMap[j - 1];
                    varMap[invVarMap[j]] = j;
                }
                invVarMap[optimalPos] = tempVar;
                varMap[invVarMap[optimalPos]] = optimalPos;
            }
		}
		return {in, total_min, total_max}; //返回DD指针和这次sifting过程中最大和最小的DD size	
	}

}
