#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include "ConjugateGradient.h"
#include "FCC.h"

int main(int argc, char **argv)
{
	SparseMatrix A;
	vector<double> b;
	SystemOpr solve;
	int N, M, num_elem;

	FILE *file[2];

	for (int f = 0; f < 2; ++f)
	{
		file[f] = fopen(argv[f+1], "w");

		if (file[f] == NULL)
		{
			printf("ERRO na abertura de arquivo!\n");
			exit(EXIT_FAILURE);
		}
	}

	solve.setType(TRIANGULAR);

	scanf("%d %d %d\n", &N, &M, &num_elem);

	A.init(N, M);
	
	for (int i = 0; i < num_elem; ++i)
	{
		double v;
		int l, c;

		scanf("%d %d %lf \n", &l, &c, &v);

		l--; c--;

		A.put(l, c, v);

		if (l != c)
			A.put(c, l, v);
	}

	A.sortRow();

	SparseMatrix I;

	I.init(A.l_size(), A.c_size());
	
	for (int k = 0; k < A.l_size(); ++k)
		I.put(k, k, 1.);


	for (int eta = -10; eta <= 10; eta += 10)
	{
		int iter = 0;
		clock_t clk;
		double sigma[] = { 0.001, 0.0001, 0.00001 }, lim = 3;

		fprintf(file[0],"eta = %d\nsigma;normal;sigma;aprox\n", eta);
		fprintf(file[1],"eta = %d\nsigma;normal;sigma;aprox\n", eta);

		for (int k = 0; k < lim; ++k)
		{
			SparseMatrix A_ = A + sigma[k] * I, G2, G;

			fprintf(file[0], "%lf;", sigma[k]);
			fprintf(file[1], "%lf;", sigma[k]);

			// teste 1: GC em A_ usando A = LL^T
			b.assign(N, 1);

			G = FCC(A, eta);
			GC_precon(G, A, b, iter, clk);

			fprintf(file[0], "%d;", iter);
			fprintf(file[1], "%lf;", (double)(clk) / CLOCKS_PER_SEC ); 
			
			printf("%d %d %lf <\n", eta, iter, (double)(clk)/ CLOCKS_PER_SEC);

			// teste 2: GC em A_ usando A_ = LL^T
			b.assign(N, 1);

			G2 = FCC(A_, eta);
			GC_precon(G2, A_, b, iter, clk);

			fprintf(file[0], "%d;", iter);
			fprintf(file[1], "%lf;", (double)(clk) / CLOCKS_PER_SEC ); 

			printf("%d %d %lf <\n", eta, iter, (double)(clk)/ CLOCKS_PER_SEC);

			// teste 3: GC usando L~ = L * D * D_^-1 
			SparseMatrix G3(A.l_size(), A.c_size());
			vector<double> D, D_;

			D_ = G2.getDiagonal();

			for (int i = 0; i < G.l_size(); ++i)
				if (G[i].r == G[i].c) 
					G3.put(G[i].r, G[i].c, G[i].v / D_[G[i].r]);	

			b.assign(N, 1);

			G3.sortRow();
			GC_precon(G3, A_, b, iter, clk);

			fprintf(file[0], "%d\n", iter);
			fprintf(file[1], "%lf\n", (double)(clk) / CLOCKS_PER_SEC ); 

			printf("%d %d %lf <\n\n", eta, iter, (double)(clk)/ CLOCKS_PER_SEC);
		}

	}
}
