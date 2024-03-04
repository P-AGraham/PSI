struct zvec
{
	MKL_INT len;
	MKL_Complex16* vec;
};

void del_zvec(struct zvec x);

#include "Level1/Level1.h"