#ifndef DDlinear_H
#define DDlinear_H


#define SWAP_MOVE 0
#define LINEAR_TRANSFORM_MOVE 1
#define INVERSE_TRANSFORM_MOVE 2

namespace dd {

    typedef struct Move	{
		short index = -1; //执行的操作的索引
		short optype; //0表示swap，1表示lt，2表示反lt
        int ddsize = -1;
	}Move;

    


}

#endif
