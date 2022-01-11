#ifndef DDlinear_H
#define DDlinear_H


#define SWAP_MOVE 0
#define LINEAR_TRANSFORM_MOVE 1
#define INVERSE_TRANSFORM_MOVE 2
#define UPDATE_MAP 3

namespace dd {

    typedef struct Move	{
		short index; //执行操作的索引
		short pos; //实际选中将要进行操作的index
		short optype; //0表示swap，1表示lt，2表示反lt
        unsigned int ddsize; //执行操作后的size
		Move() {
			this->index = -1;
			this->pos = -1;
			this->optype = -1;
			this->ddsize = 0;
		}
		bool operator==(const Move b) const  
		{ 
			return this->index == b.index && this->optype == b.optype && this->ddsize == b.ddsize && this->pos == b.pos;  
		}  
	}Move;

    


}

#endif