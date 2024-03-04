#include "Level1/Level1.c"

void del_zvec(struct zvec x)
{
	x.len = 0;
	if(x.vec!=NULL)
		mkl_free(x.vec);
	x.vec = NULL;
	return;
}