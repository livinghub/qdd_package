#include "DDpackage.h"


namespace dd {
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

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边


        short pos = -1; //选中变量的index, 电路qubit


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

            //记录初始位置
            Move curState;
            curState.ddsize = min;
            curState.index = -1;
            curState.pos = pos;
            curState.optype = -1;
            std::list<Move> curMoveV;
            curMoveV.push_back(curState);
            
            if(pos==n-1) { //选中的变量索引就是顶部
                moveDown = linearAndSiftingDown(pos, in, varMap, curMoveV);
                assert(pos == 0);
                linearAndSiftingBackward(pos, in, varMap, moveDown);                
            } else if(pos==0) {
                moveUp = linearAndSiftingUp(pos, in, varMap, curMoveV);
                assert(pos == n-1);
                linearAndSiftingBackward(pos, in, varMap, moveUp);
            } else if(pos < n/2) { // variable is in lower half -> sifting to bottom first
                moveDown = linearAndSiftingDown(pos, in, varMap, curMoveV);
                assert(pos == 0);
                moveUp = undoMoves(pos, in, varMap, moveDown);
                moveUp = linearAndSiftingUp(pos, in, varMap, moveUp);
                assert(pos == n-1);
                linearAndSiftingBackward(pos, in, varMap, moveUp);
            } else {
                moveUp = linearAndSiftingUp(pos, in, varMap, curMoveV);
                assert(pos == n-1);
                moveDown = undoMoves(pos, in, varMap, moveUp);
                moveDown = linearAndSiftingDown(pos, in, varMap, moveDown);
                assert(pos == 0);
                linearAndSiftingBackward(pos, in, varMap, moveDown);
            }

			initComputeTable();

			// there are nodes which need to renormalized
            if (unnormalizedNodes > 0) {
                // std::clog << "{" << unnormalizedNodes << "} ";
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
        
        const auto n = static_cast<short>(in.p->v + 1);
        std::list<Move> moves;
        Move move;

        moves = prevMoves;

        while (pos < n-1 ) {
            exchangeBaseCase(pos+1, in, Map); //先对pos，pos-1执行swap
            auto ex_size = size(in);
			linearInPlace(pos+1, in, Map); //再pos，pos-1执行L.T.
			auto lt_size = size(in);

            move.index = pos+1;
            move.pos = pos+1;
            
            if(ex_size <= lt_size){ //swap效果更好
				linearInPlace(pos+1, in, Map);
                move.ddsize = ex_size;
                move.optype = SWAP_MOVE;
			} else { //lt效果更好
                move.ddsize = lt_size;
                move.optype = LINEAR_TRANSFORM_MOVE;
			}
            moves.push_back(move);

            //std::clog << "↓" << ex_size << " ";
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

        // std::clog << "向下尝试";
        
        //
        std::list<Move> moves;
        Move move;

        //
        moves = prevMoves;

        while (pos > 0) {
			//std::clog << "向下尝试";
            exchangeBaseCase(pos, in, Map); //先对pos，pos-1执行swap
            auto ex_size = size(in);
			linearInPlace(pos, in, Map); //再pos，pos-1执行L.T.
			auto lt_size = size(in);

            move.index = pos;
            move.pos = pos-1;
            
            if(ex_size <= lt_size){ //swap效果更好
				//** 抵消lt
				linearInPlace(pos, in, Map);
                move.ddsize = ex_size;
                move.optype = SWAP_MOVE;


			} else { //lt效果更好
                move.ddsize = lt_size;
                move.optype = LINEAR_TRANSFORM_MOVE;
			}
            moves.push_back(move);

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

        unsigned int bSize = UINT32_MAX;
        
        for(std::list<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ++it) {
            if(it->ddsize<bSize) {
                bSize = it->ddsize;
                optimalPos = it->pos;
            }
        }  

        for(std::list<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ++it) {
            
            if(it->ddsize == bSize) {   
                optimalPos = it->pos;             
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
        for(std::list<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ++it) {            
            if(it->optype == -1) {
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
            invmove.pos = pos;
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


    //qmdd ---opSeq----> ltqmdd
    void Package::qmdd2ltqmdd(
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap, //circuit var -> DD variable
        bool Mat[MAXN][MAXN]
    ) {
        // computeMatrixProperties = Disabled;
		// Edge root{in}; //声明root边


        // static bool matInit = true;
        // if(matInit) {
        //     xorInit(varMap, Mat);
        // }
        xorInit(varMap, Mat);

        if(!opSeq.empty()) {
            
            // std::clog << std::endl << "opSequence不为空，执行转换" << std::endl;
            
            for(const short &i : opSeq) {
                // std::clog<<i<<' ';
                if(i>0) {
                    exchangeBaseCase(i, in, varMap, false);
                } else {                   
                    linearInPlace(-i, in, varMap, false); //不更新lt矩阵
                    xorLinear(-i, varMap, Mat); 
                }
            }          
        } else {
            // std::clog << std::endl << "opSequence为空" << std::endl;
        }

        
        
    }

    // ltqmdd ----(opSequence)---> qmdd
    //
    void Package::ltqmdd2qmdd(
        Edge in, 
        std::map<unsigned short, unsigned short>& varMap,
        bool Mat[MAXN][MAXN]
    ) {
        if(opSeq.empty()) return;

        std::vector<short>::reverse_iterator re_it = opSeq.rbegin();

       

        for(re_it=opSeq.rbegin(); re_it!=opSeq.rend(); ++re_it) {
            if(*re_it>0) {
                exchangeBaseCase(*re_it, in, varMap, false);
            } else {
                linearInPlace(-(*re_it), in, varMap, false); //不更新lt矩阵
                xorLinear(-(*re_it), varMap, Mat); 
            }
        }
        
        return;
    }

    Edge Package::dynTraceBack(Edge in, std::map<unsigned short, unsigned short>& varMap) {
        int mini = opSeq.size()-1;
        unsigned int miniSize = size(in);
        unsigned int curSize = miniSize;
        for(int i=mini; i>=0; --i) {
            if(opSeq[i] > 0) {
                exchangeBaseCase(opSeq[i], in);
               
            } else {
                linearInPlace(-opSeq[i], in);
                xorLinear(-opSeq[i], varMap);
            }

            curSize = size(in);
            if(curSize <= miniSize) {
                miniSize = curSize;
                mini = i;
            }
        }

        for(int i=0; i<mini; ++i) {
            if(opSeq[i] > 0) {
                exchangeBaseCase(opSeq[i], in);
               
            } else {
                linearInPlace(-opSeq[i], in);
                xorLinear(-opSeq[i], varMap);
            }
        }
        
        opSeq.resize(mini+1);

        return in;
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

		/**
        static bool matInit = true;
		if(matInit) {
            xorInit(varMap);
            matInit = false;
        }
		**/


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

