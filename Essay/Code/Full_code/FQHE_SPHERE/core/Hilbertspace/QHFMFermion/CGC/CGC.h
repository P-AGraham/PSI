

int createCG2Bdoy(struct CG2 *cg, int u);

int calCG2Body(struct CG2 *cg);

int normCG2Body_j(struct CG2 *cg, int j);

void calCG2Body_j(struct CG2 *cg, int j);

double calcg(struct CG2 *cg, int j, int x, int y, double m1, double m2);

int destoryCG2(struct CG2 *cg);

void printCG2(struct CG2 *cg);