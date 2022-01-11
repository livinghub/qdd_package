#include "DDpackage.h"


namespace dd {

    using LTorderMap = std::map<unsigned short, std::vector<unsigned short>>; // varible -> LT varible
	std::list<Move> MNULL;

    std::list<Move> findPath(std::list<Move> moves);
    std::list<Move> findPath(std::list<Move> movesA, std::list<Move> movesB);
    std::list<Move> findDirPath(std::list<Move> movesUP, std::list<Move> movesDown, short orgOps, short optOps);
    void printMoves(std::list<Move> moves);

    
    Edge Package::linearAndSiftingAux(Edge in, 
    std::map<unsigned short, unsigned short>& varMap,
    bool fg
    ) {
        //
        const auto n = static_cast<short>(in.p->v + 1); //变量个数
        //Move move;
        std::list<Move> moveUp; //list of up moves
        std::list<Move> moveDown;

		std::vector<bool> free(n, true); //记录变量是否已经被处理
		std::map<unsigned short, unsigned short> invVarMap{}; //index->var; DD qubit（变量） 到 电路 qubit（变量） 的映射
		for ( auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边


        short pos = -1; //选中变量的index

		/**/
		
		// xorInit(varMap);

		/**/


		for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;

            // std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > (unsigned short)max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos，电路变量index
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的DD变量为处理状态
            //move.ddVar = varMap[pos]; //把已经处理的DD变量记录下来
            short optimalPos = pos; 
            short originalPos = pos;

            //记录初始位置
            Move curState;
            curState.ddsize = min;
            curState.index = -1;
            curState.pos = pos;
            curState.optype = -1;
            std::list<Move> curMoveV;
            curMoveV.push_back(curState);
            
            if(pos==n-1) { //选中的变量索引就是顶部
                //
                // std::clog<<"全部向下; "<<std::endl;
                
                // std::clog<<"下: ";
                moveDown = linearAndSiftingDown(pos, in, invVarMap, curMoveV);
                // std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveDown);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
                
            } else if(pos==0) {
                //
                // std::clog<<"全部向上; "<<std::endl;
                // std::clog<<"上: ";
                moveUp = linearAndSiftingUp(pos, in, invVarMap, curMoveV);
                // std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveUp);
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
            } else if(pos < n/2) { // variable is in lower half -> sifting to bottom first
                //
                // std::clog<<"先下后上; "<<std::endl;
                std::clog<<pos<<' ';
                std::clog<<"下: ";
                moveDown = linearAndSiftingDown(pos, in, invVarMap, curMoveV);
                if(moveDown.back().index != pos+1) {
                    throw std::logic_error("检查aux-lsDown");
                }
                std::clog<<pos<<' ';
                std::clog<<"还: ";
                moveUp = undoMoves(pos, in, invVarMap, moveDown);
                if(min != size(in)) throw std::logic_error("undo 错误 ");
                std::clog<<pos<<' ';
                std::clog<<"上: ";
                moveUp = linearAndSiftingUp(pos, in, invVarMap, moveUp);
                std::clog<<pos<<' ';
                std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveUp);
                //--optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
               
            } else {
                //
                // std::clog<<"先上后下; "<<std::endl;
                
                // std::clog<<pos<<' ';
                // std::clog<<"上: ";
                moveUp = linearAndSiftingUp(pos, in, invVarMap, curMoveV);
                // if(moveUp.back().index != pos) {
                //     throw std::logic_error("检查aux-lsUp");
                // }
                
                // std::clog<<pos<<' ';
                // std::clog<<"还: ";
                moveDown = undoMoves(pos, in, invVarMap, moveUp);
                if(min != size(in)) throw std::logic_error("undo 错误 ");
                // --pos; //补偿
                
                // std::clog<<pos<<' ';
                // std::clog<<"下: ";
                moveDown = linearAndSiftingDown(pos, in, invVarMap, moveDown);
                
                // std::clog<<pos<<' ';
                // std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveDown);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
              
            }
            
            if(optimalPos == originalPos) {
                if(min != size(in)) throw std::logic_error("处理后的初始位置的size不等于原初始size");
            } else if(min < size(in)) {
                std::clog << size(in) << ' ';
                throw std::logic_error("筛选后结果更差? ");
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

			for ( auto & i : invVarMap)
			    varMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射
            
            //for debug
            // auto scheck = size(in);
            // std::clog << "    " << i << "/" << n << " size=" << scheck << " | ";
		}
		return in; //返回DD指针	
	}

    

    std::list<Move> Package::linearAndSiftingUp(
        short &pos, // variable index
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::list<Move> &prevMoves
        ) {
        //判断pos与传入的prevMoves是否冲突
        // if(!prevMoves.empty()){
        //    assert(prevMoves.back().index == pos);
        // }
        // std::clog << "向上尝试";
        
        //
        const auto n = static_cast<short>(in.p->v + 1);
        std::list<Move> moves;

        //
        moves = prevMoves;

        while (pos+1 < n ) {
			// std::clog << "向上尝试";
            exchangeBaseCase(pos+1, in, Map); //先对pos，pos-1执行swap
            auto ex_in = in;
            auto ex_size = size(in);
			linearInPlace(pos+1, in, Map); //再pos，pos-1执行L.T.
			auto lt_size = size(in);
            // std::clog << ex_size <<","<<lt_size<<" ";

            Move move;
            move.index = pos+1;
            move.pos = pos;
            move.optype = SWAP_MOVE;
            
            if(ex_size <= lt_size){ //swap效果更好
				//** 抵消lt
				// std::clog << "-";
				linearInPlace(pos+1, in, Map);
                auto lt_lt_size = size(in);
                if(lt_lt_size!=ex_size) {
                    std::clog<< lt_lt_size<<", "<<ex_size;
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
            ++pos; //变量位置上移一位 
        }

        return moves;
    }

    std::list<Move> Package::linearAndSiftingDown(
        short &pos,  
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::list<Move> &prevMoves
        ) {
        //判断pos与传入的prevMoves是否冲突
        // if(!prevMoves.empty() && prevMoves.back().index-1 != pos){
        //     std::clog<<pos<<", "<< prevMoves.back().index<<";";
        //     throw std::logic_error("传入位置错误,检查sifting down");
        // }

        // std::clog << "向下尝试";
        
        //
        std::list<Move> moves;

        //
        moves = prevMoves;

        while (pos > 0) {
			//std::clog << "向下尝试";
            exchangeBaseCase(pos, in, Map); //先对pos，pos-1执行swap
            // auto ex_in = in;
            auto ex_size = size(in);
			linearInPlace(pos, in, Map); //再pos，pos-1执行L.T.
			auto lt_size = size(in);

            Move move;
            move.index = pos;
            move.pos = pos;
            move.optype = SWAP_MOVE;
            
            if(ex_size <= lt_size){ //swap效果更好
				//** 抵消lt
				// std::clog << "-";
				linearInPlace(pos, in, Map);
                auto lt_lt_size = size(in);
                if(lt_lt_size!=ex_size) {
                    std::clog<< lt_lt_size<<", "<<ex_size;
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
        std::map<unsigned short, unsigned short>& Map, 
        std::list<Move> &moves) {
        //找到最小size
        // std::clog<<"回到最好位置 ";
        if(moves.empty()) {
            //
            throw std::logic_error("moves为空.");
        }

        unsigned int bSize = UINT32_MAX;
        
        std::list<Move>::reverse_iterator bestIt;
        for(std::list<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ++it) {
            // std::clog << "size-ind-op " << it->ddsize<<'-'<<it->index<<'-'<<it->optype<<"; ";
            if(it->ddsize<bSize) {
                bSize = it->ddsize;
                bestIt = it;
            }
        }  
        // std::clog << "best size:" << bSize<<' ' << bestIt->index<<';';
        // for( auto &i : moves) {
        //     std::clog << "size-ind-op：" << i.ddsize<<'-'<<i.index<<'-'<<i.optype<<' ';
        // }

        for(std::list<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ) {
            
            if(it->ddsize == bSize) {
                
                // std::clog<<"已是最佳位置a。";
                return 1;
            }  
            if(it->optype == LINEAR_TRANSFORM_MOVE) {
                linearInPlace(it->index, in, Map);
            }
            if(it->optype > -1) {
                exchangeBaseCase(it->index, in, Map);
            }
            if(it->optype == INVERSE_TRANSFORM_MOVE) {
                linearInPlace(it->index, in, Map);
            }
            optimalPos = it->pos;
            ++it;
            if(it->ddsize != size(in)) {
                std::clog<<size(in)<<' ';
                throw std::logic_error("还原后的size和记录上的size不一致 ");
            } 
        }

        optimalPos = -1;
        throw std::logic_error("Backward出问题! ");
        return 0;
    }

    std::list<Move> Package::undoMoves(
        short &pos, 
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::list<Move> &moves
        ) {
        // std::clog << "回到初始位置: ";

        Move invmove;
        std::list<Move> invmoves;
        unsigned int curSize = 0;
        //
        for(std::list<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ) {            
            if(it->optype == -1) {
                if(curSize != it->ddsize) throw std::logic_error("undo() 出错");
                break;
            }
            if(it->optype == SWAP_MOVE) {
                //std::clog << "ex" << it->index <<", ";
                invmove.optype = SWAP_MOVE;
                exchangeBaseCase(it->index, in, Map);
            } else if(it->optype == LINEAR_TRANSFORM_MOVE) {
                //std::clog << "lt-ex" << it->index <<", ";
                invmove.optype = INVERSE_TRANSFORM_MOVE;
                linearInPlace(it->index, in, Map);
                exchangeBaseCase(it->index, in, Map);
            } else {
                //std::clog << "ex-lt" << it->index <<", ";
                invmove.optype = LINEAR_TRANSFORM_MOVE;
                exchangeBaseCase(it->index, in, Map);
                linearInPlace(it->index, in, Map);
            }
            pos = it->pos;
            curSize = size(in);
            invmove.index = it->index;
            invmove.ddsize = curSize;
            ++it;
            if(it!=moves.rend() && curSize != it->ddsize) {
                std::clog<<"recS-curS: "<<it->ddsize<<'-'<<curSize<<' ';
                throw std::logic_error("undo后的size和记录上的size不一致 ");
            }
            
            invmoves.push_back(invmove);
        }
        return invmoves;
    }


    void printMoves(std::list<Move> moves) {
        std::cout << std::endl;
        for(auto &i : moves) {
            std::cout << "index:" << i.index << ",";
            std::cout << "size:" << i.ddsize << ",";
            std::cout << "optype:" << i.optype << "; ";
        }
        std::cout << std::endl;
    }


    //qmdd ---opSequence----> ltqmdd
    void Package::qmdd2ltqmdd(
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap, //circuit var -> DD variable
        std::list<Move> &opSeq,
        bool Mat[MAXN][MAXN]
    ) {
        std::map<unsigned short, unsigned short> invVarMap{}; //index->var; DD variable -> circuit var
		for ( auto & i : varMap)
			invVarMap[i.second] = i.first;
        computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边

        if(!opSeq.empty()) {
            
            std::clog << std::endl << "opSequence不为空，执行转换" << std::endl;
            
            for(const auto &i : opSeq) {
                if(i.optype == SWAP_MOVE) {
                    exchangeBaseCase(i.index, in, invVarMap, false);
                } else if(i.optype == LINEAR_TRANSFORM_MOVE) {                   
                    linearInPlace(i.index, in, invVarMap, false); //不更新lt矩阵
                    xorLinear(i.index, invVarMap, Mat); 
                } else { //
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
            }
            
            for ( auto & i : invVarMap)
			    varMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射
            return ;
        } else {
            std::clog << std::endl << "opSequence为空" << std::endl;
            return ;
        }
        
    }



    /***********************2022-1-6************************/


    //以交换的方式移动dd的变量
    std::list<Move> Package::exchangMoves(
    short &pos, 
    Edge in, 
    std::map<unsigned short, unsigned short>& Map, 
    std::list<Move> &moves
    ) {
        // std::clog << "回到初始位置: ";

        Move invmove;
        std::list<Move> invmoves(moves);
        unsigned int curSize = 0;
        
        //
        for(std::list<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ) {
            if(it->optype == -1) {
                break;
            }
            //std::clog << "ex" << it->index <<", ";
            
            invmove.optype = SWAP_MOVE;
            exchangeBaseCase(it->index, in, Map);
            
            pos = it->pos;
            curSize = size(in);
            // std::clog<<curSize<<' ';
            invmove.index = it->index;
            invmove.ddsize = curSize;
            ++it;
            
            
            invmoves.push_back(invmove);
        }
        return invmoves;
    }


    Edge Package::linearAndSiftingAux2(Edge in, 
    std::map<unsigned short, unsigned short>& varMap,
    bool fg
    ) {
        //
        const auto n = static_cast<short>(in.p->v + 1); //变量个数
        //Move move;
        std::list<Move> moveUp; //list of up moves
        std::list<Move> moveDown;

		std::vector<bool> free(n, true); //记录变量是否已经被处理
		std::map<unsigned short, unsigned short> invVarMap{}; //index->var; DD qubit（变量） 到 电路 qubit（变量） 的映射
		for ( auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边


        short pos = -1; //选中变量的index

		/**/
		
		// xorInit(varMap);

		/**/


		for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;

            std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > (unsigned short)max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos，电路变量index
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的DD变量为处理状态
            //move.ddVar = varMap[pos]; //把已经处理的DD变量记录下来
            short optimalPos = pos; 
            short originalPos = pos;
            
            //记录最开始状态
            Move curState;
            curState.ddsize = min;
            curState.index = pos;
            curState.optype = -1;
            std::list<Move> curMoveV;
            curMoveV.push_back(curState);
            
            if(pos==n-1) { //选中的变量索引就是顶部
                //
                // std::clog<<"全部向下; "<<std::endl;
                
                // std::clog<<"下: ";
                moveDown = linearAndSiftingDown(pos, in, invVarMap, curMoveV);
                // std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveDown);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
                
            } else if(pos==0) {
                //
                // std::clog<<"全部向上; "<<std::endl;
                // std::clog<<"上: ";
                moveUp = linearAndSiftingUp(pos, in, invVarMap, curMoveV);
                // std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveUp);
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
            } else if(pos < n/2) { // variable is in lower half -> sifting to bottom first
                //
                // std::clog<<"先下后上; "<<std::endl;
                // std::clog<<"下: ";
                moveDown = linearAndSiftingDown(pos, in, invVarMap, curMoveV);
                if(moveDown.back().index != pos+1) {
                    throw std::logic_error("检查aux-lsDown");
                }
                // std::clog<<"移: ";
                moveUp = exchangMoves(pos, in, invVarMap, moveDown);
                
                // std::clog<<"上: ";
                moveUp = linearAndSiftingUp(pos, in, invVarMap, moveUp);
                // std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveUp);
                //--optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
               
            } else {
                //
                // std::clog<<"先上后下; "<<std::endl;
                // std::clog<<"上: ";
                moveUp = linearAndSiftingUp(pos, in, invVarMap, curMoveV);
                // if(moveUp.back().index != pos) {
                //     throw std::logic_error("检查aux-lsUp");
                // }
                // std::clog<<"移: ";
                moveDown = exchangMoves(pos, in, invVarMap, moveUp);
                
                // --pos; //补偿
                // std::clog<<"下: ";
                moveDown = linearAndSiftingDown(pos, in, invVarMap, moveDown);
                // std::clog<<"好: ";
                linearAndSiftingBackward(optimalPos, in, invVarMap, moveDown);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
              
            }
            
            if(optimalPos == originalPos) {
                // if(min != size(in)) throw std::logic_error("处理后的初始位置的size不等于原初始size");
            } else if(min < size(in)) {
                std::clog << size(in) << ' ';
                throw std::logic_error("筛选后结果更差? ");
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

			for ( auto & i : invVarMap)
			    varMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射
            
            //for debug
            // auto scheck = size(in);
            // std::clog << "    " << i << "/" << n << " size=" << scheck << " | ";
		}
        
		return in; //返回DD指针	
	}
}

