#ifndef DDlinear_H
#define DDlinear_H


#define SWAP_MOVE 0
#define LINEAR_TRANSFORM_MOVE 1
#define INVERSE_TRANSFORM_MOVE 2
#define UPDATE_MAP 3

namespace dd {

    typedef struct Move	{
		//执行的动作
		short index; //执行操作的索引
		short optype; //0表示swap，1表示lt，2表示反lt
		//执行操作后的状态
		short pos; //执行操作后的位置
        unsigned int ddsize; //执行操作后的size
		Move() {
			this->index = -1;
			this->pos = -1;
			this->optype = -1;
			this->ddsize = 0;
		}
		bool operator==(const Move b) const  
		{ 
			return this->index==b.index&&this->optype==b.optype&&this->ddsize==b.ddsize&&this->pos==b.pos;  
		}  
	}Move;

	typedef struct xorNode {
		unsigned short tvar; //上层 dd变量 
		unsigned short cvar; //下层 dd变量
		xorNode(unsigned short hight_DDqubit, unsigned short low_DDqubit) {
			this->tvar = hight_DDqubit;
			this->cvar = low_DDqubit;
		}
	}xorNode;

    


}

#endif
