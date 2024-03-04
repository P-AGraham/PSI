
int createCG2Bdoy(struct CG2 *cg, int u)// SU(u1): u1 = 2 -> spin 1/2 j1 = 1/2
{
	if(u<2)
		return -1;
	if(cg->cg!=NULL)
		return -2;
	cg->u = u;
	cg->s = (double)(u-1.0)/2.0;
	cg->cg = (double **)malloc((size_t)u*sizeof(double));
	cg->head = (double *)malloc((size_t)u*u*u*sizeof(double));
	for(size_t i=0;i<(size_t)u*u*u;++i)
		cg->head[i] = -5.0;
	for(size_t i=0;i<u;++i)
		cg->cg[i] = cg->head + i*u*u; 
	return calCG2Body(cg);
}

int calCG2Body(struct CG2 *cg)
{
	if(cg->u<2)
		return -1;
	if(cg->cg==NULL||cg->head==NULL)
		return -2;
	int u = cg->u;
	for(size_t i=0;i<u;++i)
	{
		normCG2Body_j(cg, (int)i);
		calCG2Body_j(cg, (int)i);
	}
	return 0;
}

int normCG2Body_j(struct CG2 *cg, int j)
{
	double c, norm = 1.0;
	int u, m1, m2;
	double jt, mt1, mt2, mt, s;
	double *cgt = cg->cg[j];
	u = cg->u;//(u, j)
	s = cg->s;
	cgt[(u-j-1)*u+u-1] = 1.0;
	jt = (double)j;
	for(int i=1;i<u-j;++i)
	{
		m1 = u - i - 1;
		m2 = u - j - i - 1;
		mt1 = cg->s - i;
		mt2 = jt - mt1;
		mt = mt1 + mt2;
		c = -sqrt((s+mt2)*(s-mt2+1.0)/((s+mt1+1.0)*(s-mt1)));
		cgt[m2*u+m1] = cgt[(m2+1)*u+m1+1] * c;
		norm += (cgt[m2*u+m1]*cgt[m2*u+m1]);
	}
	norm = sqrt(norm);
	for(int i=0;i<u-j;++i)
	{
		m1 = u - i - 1;
		m2 = u - j - i - 1;
		cgt[m2*u+m1] /= norm;
	}
	return 0;
}

void calCG2Body_j(struct CG2 *cg, int j)
{
	int u = cg->u;
	double m1, m2, s = cg->s;
	for(int x=0;x<u;++x)
	{
		m1 = (double)x - s;
		for(int y=0;y<u;++y)
		{
			m2 = s - (double)y;
			calcg(cg, j, x, y, m1, m2);
		}
	}
	return;
}

double calcg(struct CG2 *cg, int j, int x, int y, double m1, double m2)
{
	double c1, c2, jt, s = cg->s;
	double *cgt = cg->cg[j];
	int u = cg->u;
	jt = (double)j;
	if(((fabs(m1)-s)>0.1)||((fabs(m2)-s)>0.1)||(j>u))
		return 0.0;
	if(cgt[y*u+x]>-1.1)
		return cgt[y*u+x];
	if(fabs(m1+m2)>0.1+j)
	{
		cgt[y*u+x] = 0.0;
		return 0.0;
	}
	c1 = sqrt((s+m1+1.0)*(s-m1)/((jt+m1+m2+1)*(jt-m1-m2)));
	c2 = sqrt((s+m2+1.0)*(s-m2)/((jt+m1+m2+1)*(jt-m1-m2)));
	cgt[y*u+x] = calcg(cg, j, x+1, y, m1+1, m2)*c1 + calcg(cg, j, x, y-1, m1, m2+1)*c2;
	return cgt[y*u+x];
}

int destoryCG2(struct CG2 *cg)
{
	if(cg->cg==NULL||cg->head==NULL)
		return -2;
	free(cg->cg);
	free(cg->head);
	cg->u = -1;
	return 0;
}

void printCG2(struct CG2 *cg)
{
	double *cgt;
	int u = cg->u;
	for(int j=0;j<cg->u;++j)
	{
		cgt = cg->cg[j];
		printf("************************************************************************************\n");
		printf("****************************************j=%d*****************************************\n", j);
		printf("************************************************************************************\n");
		for(int m2=0;m2<u;++m2)
		{
			for(int m1=0;m1<u;++m1)
			{
				printf("%e\t", cgt[m2*u+m1]);
			}
			printf("\n");
		}
	}
}