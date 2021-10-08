#include "DDpackage.h"


namespace dd {

	std::vector<Move> MNULL;
	//std::vector<std::vector<Move>> Movetab;

    std::vector<Move> findPath(std::vector<Move> moves);
    std::vector<Move> findPath(std::vector<Move> movesA, std::vector<Move> movesB);
    std::vector<Move> findDirPath(std::vector<Move> movesUP, std::vector<Move> movesDown, short orgOps, short optOps);
    void printMoves(std::vector<Move> moves);


    Edge Package::linearAndSiftingAux(Edge in, 
    std::map<unsigned short, unsigned short>& varMap
    ) {
        //
        const auto n = static_cast<short>(in.p->v + 1); //变量个数
        //Move move;
        std::vector<Move> moveUp; //list of up moves
        std::vector<Move> moveDown;

		std::vector<bool> free(n, true); //记录变量是否已经被处理
		std::map<unsigned short, unsigned short> invVarMap{}; //DD qubit（变量） 到 电路 qubit（变量） 的映射
		for (const auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边


        short pos = -1; //选中变量的index

		/**/
		
		xorInit(varMap);

		/**/


		for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;

            std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > (unsigned short)max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的DD变量为处理状态
            //move.ddVar = varMap[pos]; //把已经处理的DD变量记录下来
            short optimalPos = pos; 
            short originalPos = pos;
            
            if(pos==n-1) { //选中的变量索引就是顶部
                //
                // std::clog<<"全部向下; ";
                Move curState;
                curState.ddsize = min;
                curState.index = pos;
                curState.optype = -1;
                std::vector<Move> curMoveV;
                curMoveV.push_back(curState);

                moveDown = linearAndSiftingDown(pos, in, varMap, curMoveV);
                linearAndSiftingBackward(optimalPos, in, varMap, moveDown);
                --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
                
            } else if(pos==0) {
                //
                // std::clog<<"全部向上; ";
                Move curState;
                curState.ddsize = min;
                curState.index = pos;
                curState.optype = -1;
                std::vector<Move> curMoveV;
                curMoveV.push_back(curState);

                moveUp = linearAndSiftingUp(pos, in, varMap, curMoveV);
                linearAndSiftingBackward(optimalPos, in, varMap, moveUp);
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
            } else if(pos < n/2) { // variable is in lower half -> sifting to bottom first
                //
                // std::clog<<"先下后上; ";

                moveDown = linearAndSiftingDown(pos, in, varMap, MNULL);
                if(moveDown.back().index != pos+1) {
                    throw std::logic_error("检查aux-lsDown");
                }
                moveUp = undoMoves(pos, in, varMap, moveDown);
                // moveUp.push_back(curState);
                moveUp = linearAndSiftingUp(pos, in, varMap, moveUp);
                
                linearAndSiftingBackward(optimalPos, in, varMap, moveUp);
                //--optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
               
            } else {
                //
                // std::clog<<"先上后下; ";

                moveUp = linearAndSiftingUp(pos, in, varMap, MNULL);
                if(moveUp.back().index != pos) {
                    throw std::logic_error("检查aux-lsUp");
                }
                moveDown = undoMoves(pos, in, varMap, moveUp);
                --pos; //补偿
                // moveDown.push_back(curState);
                moveDown = linearAndSiftingDown(pos, in, varMap, moveDown);
  
                linearAndSiftingBackward(optimalPos, in, varMap, moveDown);
                --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
              
            }
            
            if(optimalPos == originalPos) {
                if(min != size(in)) throw std::logic_error("处理后的初始位置的size不等于原初始size");
            }
                    

            if(0) {
                    std::clog<<std::endl;
                    std::clog<<"--------------------";
                    // printMoves(moveUp);
                    // printMoves(moveDown);
                    // printMoves(findDirPath(moveUp, moveDown, originalPos, optimalPos));
                    printLTPath(LTpath);
                    std::clog<<"--------------------";
                    std::clog<<std::endl;
            }
            

            assert(optimalPos>=0);

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

            //for debug
            // auto scheck = size(in);
            // std::clog << "    " << i << "/" << n << " size=" << scheck << " | ";
		}
		return in; //返回DD指针	
	}

    

    std::vector<Move> Package::linearAndSiftingUp(
        short &pos, 
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap, 
        std::vector<Move> prevMoves
        ) {
        //判断pos与传入的prevMoves是否冲突
        // if(!prevMoves.empty()){
        //    assert(prevMoves.back().index == pos);
        // }
        // std::clog << "向上尝试";
        
        //
        const auto n = static_cast<short>(in.p->v + 1);
        std::vector<Move> moves;

        //
        moves = prevMoves;

        while (pos < n - 1) {
			// std::clog << "向上尝试";
            exchangeBaseCase(pos+1, in); //先对pos，pos-1执行swap
            auto ex_in = in;
            auto ex_size = size(in);
			linearInPlace(pos+1, in, varMap); //再pos，pos-1执行L.T.
			auto lt_size = size(in);

            Move move;
            move.index = pos+1;
            move.optype = SWAP_MOVE;
            
            if(ex_size <= lt_size){ //swap效果更好
				//** 抵消lt
				// std::clog << "-";
				linearInPlace(pos+1, in, varMap);
                auto lt_lt_size = size(in);
                if(lt_lt_size!=ex_size || !equals(ex_in,in)) {
                    throw std::logic_error("lt-lt 后,非还原状态");
                }
                move.ddsize = ex_size;

                //debug
                // std::clog<<"swap better"<<std::endl;
			} else { //lt效果更好
                //debug
                // std::clog<<"lt better"<<std::endl;

                move.ddsize = lt_size;
                move.optype = LINEAR_TRANSFORM_MOVE;
			}
            moves.push_back(move);

            //std::clog << "↓" << ex_size << " ";
            assert(is_locally_consistent_dd(in));
            ++pos; //变量位置下移一位 
        }

        return moves;
    }

    std::vector<Move> Package::linearAndSiftingDown(
        short &pos,  
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap, 
        std::vector<Move> prevMoves
        ) {
        //判断pos与传入的prevMoves是否冲突
        // if(!prevMoves.empty() && prevMoves.back().index-1 != pos){
        //     std::clog<<pos<<", "<< prevMoves.back().index<<";";
        //     throw std::logic_error("传入位置错误,检查sifting down");
        // }

        // std::clog << "向下尝试";
        
        //
        std::vector<Move> moves;

        //
        moves = prevMoves;

        while (pos > 0) {
			//std::clog << "向下尝试";
            exchangeBaseCase(pos, in); //先对pos，pos-1执行swap
            auto ex_in = in;
            auto ex_size = size(in);
			linearInPlace(pos, in, varMap); //再pos，pos-1执行L.T.
			auto lt_size = size(in);

            Move move;
            move.index = pos;
            move.optype = SWAP_MOVE;
            
            if(ex_size <= lt_size){ //swap效果更好
				//** 抵消lt
				// std::clog << "-";
				linearInPlace(pos, in, varMap);
                auto lt_lt_size = size(in);
                if(lt_lt_size!=ex_size || !equals(ex_in,in)) {
                    throw std::logic_error("lt-lt 后,非还原状态");
                }
                move.ddsize = ex_size;

                //debug
                // std::clog<<"swap better"<<std::endl;
			} else { //lt效果更好
                move.ddsize = lt_size;
                move.optype = LINEAR_TRANSFORM_MOVE;

                //debug
                // std::clog<<"lt better"<<std::endl;
			}
            moves.push_back(move);

            assert(is_locally_consistent_dd(in));
            --pos; //变量位置下移一位 
        }

        return moves;

    }

    int Package::linearAndSiftingBackward(
        short &optimalPos, 
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap, 
        std::vector<Move> moves) {
        //找到最小size
        // std::clog<<"回到最好位置 ";
        if(moves.empty()) {
            //
            throw std::logic_error("moves为空！");
        }

        int size=moves.front().ddsize;
        for(std::vector<Move>::iterator it=moves.begin(); it!=moves.end(); ++it) {
            if(it->ddsize < size) {
                size = it->ddsize;
            }
        }
        // std::clog << "best size：" << size << ", ";

        for(std::vector<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ++it) {
            if(it->ddsize == size) {
                optimalPos = it->index;
                // std::clog<<"已是最佳位置。";
                return 1;
            }
            if(it->optype == LINEAR_TRANSFORM_MOVE) {
                linearInPlace(it->index, in, varMap);
            }
            if(it->optype > -1) {
                exchangeBaseCase(it->index, in);
            }
            if(it->optype == INVERSE_TRANSFORM_MOVE) {
                linearInPlace(it->index, in, varMap);
            }
        }

        optimalPos = -1;
        return 0;
    }

    std::vector<Move> Package::undoMoves(
        short &pos, 
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap, 
        std::vector<Move> moves
        ) {
        // std::clog << "回到初始位置: ";

        Move invmove;
        std::vector<Move> invmoves;
        //
        for(std::vector<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ++it) {
            pos = it->index;
            if(it->optype == -1) break;
            if(it->optype == SWAP_MOVE) {
                //std::clog << "ex" << it->index <<", ";
                invmove.optype = SWAP_MOVE;
                exchangeBaseCase(it->index, in);
            } else if(it->optype == LINEAR_TRANSFORM_MOVE) {
                //std::clog << "lt-ex" << it->index <<", ";
                invmove.optype = INVERSE_TRANSFORM_MOVE;
                linearInPlace(it->index, in, varMap);
                exchangeBaseCase(it->index, in);
            } else {
                //std::clog << "ex-lt" << it->index <<", ";
                invmove.optype = LINEAR_TRANSFORM_MOVE;
                exchangeBaseCase(it->index, in);
                linearInPlace(it->index, in, varMap);
            }
            invmove.index = it->index;
            invmove.ddsize = size(in);
            invmoves.push_back(invmove);
        }
        return invmoves;
    }

    void Package::printLTPath(std::vector<std::pair<short, short>> LTpath) {
        std::clog<<std::endl;
        for(auto &i : LTpath) {
            std::clog << i.first << "-" << i.second << ";    ";
        }
        std::clog<<std::endl;
    }


    // undo
    bool Package::findLTPath(
        std::vector<Move> movesUP, 
        std::vector<Move> movesDown, 
        std::map<unsigned short, unsigned short>& varMap,
        short orgOps, 
        short optOps
        ) {
        //
        
        std::vector<Move>::iterator begenit, endit;
        std::pair<short, short> LTentry;
        

        //找到最好的位置
        if(optOps==orgOps) {
            //
            return 0;
        } else if(optOps>orgOps) {
            //
            for(begenit=movesUP.begin(); begenit->index<orgOps; ++begenit);
            for(endit=begenit; endit->index<=optOps;++endit) {
                if(endit->optype == LINEAR_TRANSFORM_MOVE) {
                    // LTentry = std::make_pair(endit->index-1, endit->index);
                    // LTpath.push_back(LTentry);
                    LTentry = std::make_pair(varMap[endit->index], varMap[endit->index-1]);
                    LTpath.push_back(LTentry);
                }
            }
            return 1;
        } else {
            //
            for(begenit=movesDown.begin(); begenit->index>orgOps; ++begenit);
            for(endit=begenit; endit->index>=optOps+1; ++endit) {
                if(endit->optype == LINEAR_TRANSFORM_MOVE) {
                    // LTentry = std::make_pair(endit->index, endit->index-1);
                    // LTpath.push_back(LTentry);
                    LTentry = std::make_pair(varMap[endit->index], varMap[endit->index-1]);
                    LTpath.push_back(LTentry);
                }
            }
            return 1;
        }
        
    }

    //输入向上移动和向下移动两段移动记录,找到由原始位置到最好位置的一段记录
    std::vector<Move> findDirPath(std::vector<Move> movesUP, std::vector<Move> movesDown, short orgOps, short optOps) {
        //
        
        std::vector<Move>::iterator begenit, endit;

        //找到最好的位置
        if(optOps==orgOps) {
            //
            std::vector<Move> movePath;
            return movePath;
        } else if(optOps>orgOps) {
            //
            for(endit=movesUP.end(); endit->index!=optOps; --endit);
                //
            for(begenit=endit; begenit->index!=orgOps+1; --begenit);
            std::vector<Move> movePath(begenit,++endit);
            return movePath;
        } else {
            //
            for(begenit=movesDown.begin(); begenit->index!=orgOps; ++begenit);
            for(endit=begenit; endit->index!=optOps+1; ++endit);
            std::vector<Move> movePath(begenit,++endit);
            return movePath;
        }
        
    }

    //输入向某个方向尝试的移动记录,找到由原始位置到最好位置的一段记录
    std::vector<Move> findPath(std::vector<Move> moves) {
        //
        
        std::vector<Move>::iterator it;
        int minSize = moves.front().ddsize;

        //找到最好的位置
        for(it=moves.begin(); it!=moves.end(); ++it) {
            if(it->ddsize < minSize) {
                minSize = it->ddsize;
            }
        }

        //连同操作顺序一起还原
        for(it=moves.end()-1; it->ddsize!=minSize; --it);

        std::vector<Move> movePath(moves.begin(), ++it);
        return movePath;
    }

    //输入向上尝试和向下尝试的两段移动记录,找到由原始位置到最好位置的一段记录
    std::vector<Move> findPath(std::vector<Move> movesA, std::vector<Move> movesB) {
        //
        int minA = movesA.front().ddsize, minB = movesB.back().ddsize;
        std::vector<Move>::iterator minitA, minitB;

        for(minitA = movesA.begin(); minitA!=movesA.end(); ++minitA) {
            if(minitA->ddsize < minA) {
                minA = minitA->ddsize;
            }
        }

        for(minitB = movesB.begin()+movesA.size(); minitB != movesB.end(); ++minitB) {
            if(minitB->ddsize < minB) {
                minB = minitB->ddsize;
            }
        }

        if(minA <= minB) {
            for(minitA = movesA.begin(); minitA->ddsize != minA; ++minitA);
            
            std::vector<Move> movePath(movesA.begin(), ++minitA);
            return movePath;
        } else {
            for(minitB = movesB.end()-1; minitB->ddsize != minB; --minitB);
           
            std::vector<Move> movePath(movesB.begin()+movesA.size(), ++minitB);
            return movePath;
        }
    }

    void printMoves(std::vector<Move> moves) {
        std::cout << std::endl;
        for(auto &i : moves) {
            std::cout << "index:" << i.index << ",";
            std::cout << "size:" << i.ddsize << ",";
            std::cout << "optype:" << i.optype << "; ";
        }
        std::cout << std::endl;
    }

    // 打印movetab
    void Package::printMovetab(std::vector<std::vector<Move>> Movetab) {
        for(auto i:Movetab) {
            for(auto j:i) {
                std::cout<<"索引:"<<j.index<<",操作:"<<j.optype<<" ;";
            }
        }
        std::cout<<std::endl;
    }

    //由qmdd to ltqmdd
	//输入qmdd and 变量映射表 Movetab；输出ltqmdd
    Edge Package::qmdd2ltqmdd(
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap,
        std::vector<std::vector<Move>> Movetab
    ) {
        std::map<unsigned short, unsigned short> invVarMap{}; //DD qubit（变量） 到 电路 qubit（变量） 的映射
		for (const auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

        if(!Movetab.empty()) {
            std::clog << std::endl << "movetab不为空，执行转换" << std::endl;
            for (auto &moves : Movetab) {
                //
                for(auto &i : moves) {
                    if(i.optype == SWAP_MOVE) {
                        exchangeBaseCase(i.index, in);
                        auto tempVar = invVarMap[i.index];
                        invVarMap[i.index] = invVarMap[i.index-1];
                        invVarMap[i.index-1] = tempVar;

                    } else {
                        exchangeBaseCase(i.index, in);
                        linearInPlace(i.index, in, varMap, false); //搞多一个矩阵
                        auto tempVar = invVarMap[i.index];
                        invVarMap[i.index] = invVarMap[i.index-1];
                        invVarMap[i.index-1] = tempVar;
                    }
                }
                std::clog << ";    ";
            }
            for (const auto & i : invVarMap)
			    varMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射
            return in;
        } else {
            std::clog << std::endl << "movetab为空" << std::endl;
            return in;
        }
        
    }


// 根据先前执行的linear sifting记录下的 dd变量的lt操作, 对新的qmdd -> ltqmdd
// 根据LTpath
 Edge Package::qmdd2ltqmdd(
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap
    ) {
        std::map<unsigned short, unsigned short> invVarMap{};
        computeMatrixProperties = Disabled;
        Edge root{in}; //声明root边
        for (const auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射
        
        if(!LTpath.empty()) {
            for(auto &i : LTpath) {
                if(invVarMap[i.first] > invVarMap[i.second]) {
                    while(invVarMap[i.first] != invVarMap[i.second]+1) {
                        exchangeBaseCase(invVarMap[i.first], in);
                        auto temp = i.first;
                        varMap[invVarMap[i.first]] = varMap[invVarMap[i.first]-1];
                        varMap[invVarMap[i.first]-1] = temp;
                        --invVarMap[i.first];
                        
                    }

                    linearInPlace(invVarMap[i.first], in);

                    // for(const auto& Q: varMap) {
                    //     std::clog <<"\t" << Q.first << ": " << Q.second << std::endl;
                    // }
                } else {
                    while (invVarMap[i.second] != invVarMap[i.first])
                    {
                        exchangeBaseCase(invVarMap[i.second], in);
                        auto temp = i.second;
                        varMap[invVarMap[i.second]] = varMap[invVarMap[i.second]-1];
                        varMap[invVarMap[i.second]-1] = temp;
                        --invVarMap[i.second];
                        
                    }
                    linearInPlace(invVarMap[i.first]+1, in);

                    //debug
                    // for(const auto& Q: varMap) {
                    //     std::clog <<"\t" << Q.first << ": " << Q.second << std::endl;
                    // }
                }
                for (const auto & j : varMap)
			                invVarMap[j.second] = j.first; //DD qubit 到 电路 qubit 的映射

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
            }
            
            return in;
        } else {
            return in;
        }
        // for (const auto & i : invVarMap)
		// 	    varMap[i.second] = i.first; //cir qubit 到 dd qubit 的映射

    }

    //获取当前ltqmdd的movetab
    //其实movetab是公有成员,可以不用这个函数
    std::vector<std::vector<Move>> Package::getMovetab() {
        return Movetab;
    }

    // ltqmdd ----(movetab)---> qmdd
    //
    Edge Package::ltqmdd2qmdd(
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap,
        std::vector<std::vector<Move>> Movetab
    ) {
        assert(!Movetab.empty());

        std::vector<std::vector<Move>>::reverse_iterator re_it = Movetab.rbegin();
        std::vector<Move>::reverse_iterator re_it2;

        std::map<unsigned short, unsigned short> invVarMap{}; //DD qubit（变量） 到 电路 qubit（变量） 的映射
		for (const auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

        for(; re_it!=Movetab.rend(); ++re_it) {
            for(re_it2 = re_it->rbegin(); re_it2!=re_it->rend(); ++re_it2) {
                if(re_it2->optype == SWAP_MOVE) {
                        exchangeBaseCase(re_it2->index, in);
                        auto tempVar = invVarMap[re_it2->index];
                        invVarMap[re_it2->index] = invVarMap[re_it2->index-1];
                        invVarMap[re_it2->index-1] = tempVar;

                } else {
                    //exchangeBaseCase(re_it2->index, in);
                    linearInPlace(re_it2->index, in, varMap, false); //不刷新矩阵
                    exchangeBaseCase(re_it2->index, in);
                    auto tempVar = invVarMap[re_it2->index];
                    invVarMap[re_it2->index] = invVarMap[re_it2->index-1];
                    invVarMap[re_it2->index-1] = tempVar;
                }
            }
        }
        for (const auto & i : invVarMap)
			    varMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射
        return in;
    }

}

