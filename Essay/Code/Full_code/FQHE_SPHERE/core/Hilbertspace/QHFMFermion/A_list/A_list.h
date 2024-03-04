
void cal2BodyAlist(struct CG2 *cg, double *pseudo, double *A_list2, int num_thread);

void setCoulomb(double *pseudo, size_t N_phi);

// O00 = \sum_l \sum_{m1,m2,m} o_l T^0_0(l)
void calO00Alist(struct CG2 *cg, double *olist, double *A_list2, int num_thread);

void cal_nlm_nlm_Alist(double *A_list2, int l, int m, int N_o, int num_thread);

void cal_nlm_nlm_Alist(double *A_list2, double *olist, int l, int N_o, int num_thread);

void cal_nl1m_nl2m_Alist(double *A_list2, int l_1, int l_2, int m, int N_o, int num_thread);

void cal_delta_Alist(double *A_list2, int l, int m, int N_o, int num_thread);

void cal_deltadelta_lalb_Alist(double *A_list2, int l, int m, int la, int lb, int N_o, int num_thread);

void cal_nlm_Alist(double *A_list2, int l, int m, int N_o, int num_thread);