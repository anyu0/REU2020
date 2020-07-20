//---------------------------------------------------------------------
// An exact sampler for the qRacah ensemble
//---------------------------------------------------------------------
public class SimulateQR {

    public double q, k;          // the q and kappa parameters
    public int N, T, S, M;          // the sizes of the hexagon
    public double[] qM;          // matrix of q^m
    public double[] logQM;          // matrix of m*log(q)
    public double[] logQMm1;          // matrix of logarithms of 1 - q^m
    public double[] logQMm1k;          // matrix of logarithms of 1 - k*q^m

//*************************************************************************
// Basic functions for the class: Constructor, sampler, and helper
// and helper functions for both
//*************************************************************************

    //---------------------------------------------------------------------
    // Initialize the parameters and make a matrix qM
    //---------------------------------------------------------------------
    public SimulateQR(double q, double k, int N, int T, int S) {
        this.N = N;
        this.S = S;
        this.T = T;
        this.q = q;
        this.k = k;
        this.M = 5000; // a large constant
        this.qM = new double[4 * M + 1];
        this.logQM = new double[4 * M + 1];
        this.logQMm1 = new double[2 * M + 1];
        this.logQMm1k = new double[4 * M + 1];
        qM[2 * M] = 1;
        logQMm1k[2 * M] = Math.log(1 - k);
        for (int i = 1; i < 2 * M; i++) {
            logQM[i + 2 * M] = logQM[i - 1 + 2 * M] + Math.log(q);
            logQM[-i + 2 * M] = logQM[-i + 1 + 2 * M] - Math.log(q);

            qM[i + 2 * M] = qM[i - 1 + 2 * M] * q;
            qM[-i + 2 * M] = qM[-i + 1 + 2 * M] / q;

            logQMm1[i] = Math.log(1 - qM[i + 2 * M]);
            logQMm1k[2 * M + i] = Math.log(1 - k * qM[i + 2 * M]);
            logQMm1k[2 * M - i] = Math.log(1 - k * qM[-i + 2 * M]);
        }


    }

    //---------------------------------------------------------------------
    // A helper function qM - gives powers of q
    //---------------------------------------------------------------------
    private double qM(int x) {
        return qM[x + 2 * M];
    }

    //---------------------------------------------------------------------
    // A helper function logQM - gives logs of q
    //---------------------------------------------------------------------
    private double logQM(int x) {
        return logQM[x + 2 * M];
    }

    //---------------------------------------------------------------------
    // A helper function logQM - gives logs of q
    //---------------------------------------------------------------------
    private double logQMm1(int x) {
        return logQMm1[x];
    }

    //---------------------------------------------------------------------
    // A helper function logQM - gives logs of q
    //---------------------------------------------------------------------
    private double logQMm1k(int x) {
        return logQMm1k[x + 2 * M];
    }

    //---------------------------------------------------------------------
    // A helper function p - gives the "little p" from Section 6.2 in BGR
    //---------------------------------------------------------------------
    private double p(int x, int t, int s) {
        if ((x + T - t - s - 1) == 0) return 0;

        double p = (1 - qM(x + T - t - s - 1)) * (1 - k * qM(x - s - t - 1));

        p = p * (1 - k * qM(2 * x - t - s + 1)) / (qM(T - t - s - 1) * (1 - qM(x + 1)));
        p = p / ((1 - k * qM(x - T + 1)) * (1 - k * qM(2 * x - t - s - 1)));


        return p;
    }

    //---------------------------------------------------------------------
    // A helper function p - gives the "little p" from Section 6.2 in BGR
    //---------------------------------------------------------------------
    private double logP(int x, int t, int s) {

        double p = logQMm1(x + T - t - s - 1) + logQMm1k(x - s - t - 1);

        p = p + logQMm1k(2 * x - t - s + 1) - logQM(T - t - s - 1) - logQMm1(x + 1);
        p = p - logQMm1k(x - T + 1) - logQMm1k(2 * x - t - s - 1);
        return p;
    }

    //---------------------------------------------------------------------
    // A helper function D - gives the distribution from eq (14) in BGR
    //---------------------------------------------------------------------
    public double[] D(int x, int t, int s, int n) {
        double[] D = new double[n + 1];

        D[0] = 0;
        int stop = 1;
        boolean check = true;
        for (int i = 1; i < (n + 1); i++) {
            if (((x + i - 1) + T - t - (s - 1) - 1) == 0) {
                check = false;
                stop = i;
            }

            if (check) {
                D[i] = D[i - 1] + logP(x + i - 1, t, s - 1);
            }
        }
        if (check) {
            stop = n + 1;
        }
        // center D
        double max = 0;
        for (int i = 0; i < stop; i++) {
            if (D[i] > max) {
                max = D[i];
            }
        }

        double sum = 0;
        for (int i = 0; i < stop; i++) {
            if ((D[i] - max) < -50) {
                D[i] = 0;
            } else {
                D[i] = Math.exp(D[i] - max);
            }
            sum = sum + D[i];
        }
        for (int i = 0; i < stop; i++) {
            D[i] = D[i] / sum;
        }

//        D[0] = 1;
//        int stop = 1;
//        boolean check = true;
//        for (int i = 1; i < (n+1); i++){
//            if (((x+i-1) + T - t - (s-1) - 1) == 0){
//                check = false;
//                stop = i;
//            }
//
//            if (check){
//                D[i] = D[i-1]*p(x+i-1, t,s-1);;
//            }
//        }
//
//        if (check){ stop = n+1;}
//        // center D
//
//
//        double sum = 0;
//        for (int i = 0; i < stop; i++){
//            sum = sum + D[i];
//        }
//


        return D;

    }

    //---------------------------------------------------------------------
    // Simulate a new tiling of an N, S, T hexagon with the q,k parameters
    //---------------------------------------------------------------------
    public Tiling sample() {

        Tiling tile = new Tiling(N, T, 0);
        // We increase s until we reach S
        for (int s = 1; s < (S + 1); s++) {
            for (int t = 0; t < T; t++) {
                int[] x = new int[N];    // X(t)
                int[] y = new int[N];    // Y(t-1)

                int[] blocks = new int[N]; // stores the sizes of blocks
                int[] index = new int[N]; // stores smallest number in blocks

                int size = 0;           // dummy variables
                int count = 0;

                for (int i = 0; i < N; i++) {
                    x[i] = tile.particles[i][t + 1];
                    y[i] = tile.particles[i][t];
                }

                for (int i = 0; i < N; i++) {
                    if (x[i] == y[i]) {
                        size = 1;
                        for (int j = 0; j < (N - i - 1); j++) {
                            if (x[i + j + 1] != x[i + j] + 1) break;
                            if (y[i + j + 1] != y[i + j] + 1) break;
                            size++;
                        }
                        index[count] = i;
                        blocks[count] = size;
                        i = i + size - 1;
                        count++;
                    } else if ((x[i] - y[i]) == 1) {
                        tile.particles[i][t + 1] = x[i];
                    } else if ((x[i] - y[i]) == -1) {
                        tile.particles[i][t + 1] = y[i];
                    }

                }

                // Sample each of the blocks and update z
                for (int i = 0; i < count; i++) {
                    double[] d = D(x[index[i]], t, s, blocks[i]);

                    double u = Math.random();
                    double sum = 0;
                    for (int j = 0; j < blocks[i]; j++) {
                        sum = sum + d[j];
                        if (u < sum) {
                            tile.particles[index[i] + j][t + 1] = x[index[i] + j] + 1;
                        } else {
                            tile.particles[index[i] + j][t + 1] = x[index[i] + j];
                        }
                    }
                }
            }

            tile.S = s;
            tile.c = s;
            tile.b = tile.b - 1;
        }
        return tile;
    }

//*************************************************************************
// Specialized functions for the class: tester of sampler
//*************************************************************************


    //---------------------------------------------------------------------
    // Test whether the sampling algorithm is correct
    // idea is to sample many tilings of N = 2, S = 1, T = 3 and see if
    // the percentage of each type of tiling match the probabilities
    //---------------------------------------------------------------------
    public void testSample() {

        // the weights of a horizontal lozenge
        double[] w = new double[3];
        w[0] = k / q - q;
        w[1] = k - 1;
        w[2] = k * q - 1 / q;

        // the probabilities of the 6 possible configurations
        double[] prob = new double[6];
        prob[0] = w[0] * w[0];
        prob[1] = w[1] * w[0];
        prob[2] = w[1] * w[1];
        prob[3] = w[2] * w[0];
        prob[4] = w[2] * w[1];
        prob[5] = w[2] * w[2];

        double sum = 0;
        for (int i = 0; i < 6; i++) {
            sum = sum + prob[i];
        }
        for (int i = 0; i < 6; i++) {
            prob[i] = prob[i] / sum;
        }

        // sample iter many tilings and record their frequencies
        int iter = 1000000;
        int[] freq = new int[6];
        for (int i = 0; i < iter; i++) {
            Tiling tile = this.sample();
            int type = 0;  // stores the type of the tiling - 6 possibilities
            if (tile.particles[1][1] == 2) {
                if (tile.particles[0][1] == 0) {
                    if (tile.particles[0][2] == 0) {
                        type = 2;
                    } else {
                        type = 1;
                    }
                } else {
                    type = 0;
                }
            } else {
                if (tile.particles[1][2] == 2) {
                    if (tile.particles[0][2] == 0) {
                        type = 4;
                    } else {
                        type = 3;
                    }
                } else {
                    type = 5;
                }
            }
            freq[type]++;
        }
        for (int i = 0; i < 6; i++) {
            StdOut.println(prob[i] + " " + freq[i] * 1.0 / iter);
        }
    }


    //---------------------------------------------------------------------
    // Test whether the sampling algorithm is correct
    // idea is to sample many tilings of N = 4, S = 1, T = 2 and see if
    // the percentage of each type of tiling match the probabilities
    //---------------------------------------------------------------------
    public void testSample2() {

        // the weights of a horizontal lozenge
        double[] w = new double[5];
        w[0] = k / q - q;
        w[1] = k - 1;
        w[2] = k * q - 1 / q;
        w[3] = k * q * q - 1 / (q * q);
        w[4] = k * q * q * q - 1 / (q * q * q);

        // the probabilities of the 4 possible configurations
        double[] prob = new double[5];
        prob[0] = w[0];
        prob[1] = w[1];
        prob[2] = w[2];
        prob[3] = w[3];
        prob[4] = w[4];

        double sum = 0;
        for (int i = 0; i < 5; i++) {
            sum = sum + prob[i];
        }
        for (int i = 0; i < 5; i++) {
            prob[i] = prob[i] / sum;
        }

        // sample iter many tilings and record their frequencies
        int iter = 1000000;
        int[] freq = new int[5];
        for (int i = 0; i < iter; i++) {
            Tiling tile = this.sample();
            int type = 0;  // stores the type of the tiling - 4 possibilities
            for (int j = 0; j < 4; j++) {
                if (tile.particles[j][1] == j) {
                    type++;
                }
            }

            freq[type]++;
        }
        for (int i = 0; i < 5; i++) {
            StdOut.println(prob[i] + " " + freq[i] * 1.0 / iter);
        }
    }

    //---------------------------------------------------------------------
    // Test whether the sampling algorithm is correct
    // idea is to sample many tilings of N = 2, S = 2, T = 3 and see if
    // the percentage of each type of tiling match the probabilities
    //---------------------------------------------------------------------
    public void testSample3() {

        // the weights of a horizontal lozenge
        double[] w = new double[3];
        w[0] = k / q - q;
        w[1] = k - 1;
        w[2] = k * q - 1 / q;

        // the probabilities of the 6 possible configurations
        double[] prob = new double[6];
        prob[0] = w[2] * w[2];
        prob[1] = w[1] * w[2];
        prob[2] = w[2] * w[0];
        prob[3] = w[1] * w[1];
        prob[4] = w[1] * w[0];
        prob[5] = w[0] * w[0];

        double sum = 0;
        for (int i = 0; i < 6; i++) {
            sum = sum + prob[i];
        }
        for (int i = 0; i < 6; i++) {
            prob[i] = prob[i] / sum;
        }

        // sample iter many tilings and record their frequencies
        int iter = 1000000;
        int[] freq = new int[6];
        int[] freq2 = new int[3];
        for (int i = 0; i < iter; i++) {
            Tiling tile = this.sample();
            int type = 0;  // stores the type of the tiling - 6 possibilities
            if (tile.particles[0][1] == 0) {
                if (tile.particles[1][1] == 2) {
                    if (tile.particles[1][2] == 2) {
                        type = 1;
                    } else {
                        type = 3;
                    }
                } else {
                    type = 0;
                }
            } else {
                if (tile.particles[0][2] == 1) {
                    if (tile.particles[1][2] == 3) {
                        type = 4;
                    } else {
                        type = 2;
                    }
                } else {
                    type = 5;
                }
            }
            freq[type]++;

            type = 0;
            if (tile.particles[1][1] == 2) {
                if (tile.particles[0][1] == 0) {
                    type = 1;
                } else {
                    type = 2;
                }
            } else {
                type = 0;
            }
            freq2[type]++;
        }
        for (int i = 0; i < 6; i++) {
            StdOut.println(prob[i] + " " + freq[i] * 1.0 / iter);
        }
        for (int i = 0; i < 3; i++) {
            StdOut.println(freq2[i] * 1.0 / iter);
        }
    }


//*************************************************************************
// Specialized functions for the paper: limit shape
//*************************************************************************

    //---------------------------------------------------------------------
    // Sample limit shape vs actual limit shape
    //---------------------------------------------------------------------
    public void testLS() {
        int M = Math.max(N + S, T);
        StdDraw.setXscale(1, qM(-M) * (1 + k));
        StdDraw.setYscale(1, qM(-M) * (1 + k));

        // sample iter many tilings and the top and bottom hole position
        int iter = 1;
        double[][][] pts = new double[iter][T + 1][4];
        for (int i = 0; i < iter; i++) {
            Tiling tile = this.sample();
            int[][] coord = tile.boundary();
            for (int j = 1; j < T; j++) {
                double x = qM(coord[j][0]);
                double y = qM(coord[j][1]);

                pts[i][j][0] = 1 / x + x * k * qM(-S) * qM(-j);
                pts[i][j][1] = 1 / y + y * k * qM(-S) * qM(-j);
            }
        }

        double[] avx = new double[T + 1];
        double[] avy = new double[T + 1];
        for (int i = 0; i < iter; i++) {
            for (int j = 1; j < T; j++) {
                avx[j] = avx[j] + pts[i][j][0];
                avy[j] = avy[j] + pts[i][j][1];
            }
        }

        // draw empirical limit shape
        for (int j = 1; j < T; j++) {
            avx[j] = avx[j] / iter;
            avy[j] = avy[j] / iter;
            StdDraw.filledCircle(qM(-j), avx[j], qM(-M) * 0.003);
            StdDraw.filledCircle(qM(-j), avy[j], qM(-M) * 0.003);
        }

        // draw theoretical limit shape
        StdDraw.setPenColor(StdDraw.BOOK_RED);
        double T1 = qM(-T);
        double S1 = qM(-S);
        double q1 = qM(N);

        for (int j = 1; j < T; j++) {
            double x1 = qM(-j);
            double D = q1 * Math.sqrt(16 * (1 - q1) * (1 - k * q1) * (x1 - 1) * (T1 - q1) * (1 - k * T1 * q1) * (T1 - x1) * (S1 - 1) * (T1 - S1));
            double a1 = q1 * q1 * (T1 - 1) * (T1 - 1);
            double b1 = 2 * (T1 - 1) * (S1 * x1 - T1 + (2 - S1 - x1) * q1);
            b1 = b1 - 4 * (T1 - q1) * (S1 - 1) * (x1 - 1);
            b1 = b1 * q1 - 2 * q1 * k * (S1 * T1 * x1 * q1 * q1 + S1 * T1 * T1 * q1 - 2 * S1 * T1 * x1 * q1 - 2 * S1 * T1 * q1 * q1);
            b1 = b1 - 2 * q1 * k * (S1 * x1 * q1 * q1 + T1 * T1 * x1 * q1 + T1 * T1 * q1 * q1 - 2 * T1 * x1 * q1 * q1);
            b1 = b1 - 2 * q1 * k * (S1 * T1 * q1 - 2 * T1 * T1 * q1 + T1 * x1 * q1 + T1 * q1 * q1);

            double top = (-b1 + D) / (2 * a1);
            double bot = (-b1 - D) / (2 * a1);

            StdDraw.filledCircle(qM(-j), top, qM(-M) * 0.003);
            StdDraw.filledCircle(qM(-j), bot, qM(-M) * 0.003);
        }


    }

    //---------------------------------------------------------------------
    // Main method: put here if you want to do something specific
    //---------------------------------------------------------------------
    public static void main(String[] args) {
//        double q = 0.995;
//        double k = 0.1;
//        int N = 200;
//        int T = 400;
//        int S = 200;
//        int iter = 500;
//        StdOut.println(q + " " + k + " " +N + " " + T + " " + S);
//        SimulateQR qr = new SimulateQR(q,k,N,T,S);
//
//        for (int i = 0; i < iter; i++){
//            Tiling tile = qr.sample();
//            for (int j = 0; j < N; j++){
//                for (int t= 0; t < (T+1); t++){
//                    StdOut.print(tile.particles[j][t] + " ");
//                }
//                StdOut.println();
//            }
//            StdOut.println();
//        }

        double q = 0.999999;
        double k = 0;
        int N = 30;
        int T = 50;
        int S = 20;
        SimulateQR qr = new SimulateQR(q, k, N, T, S);
        Tiling tile = qr.sample();
        tile.drawLines();
    }
}
