//---------------------------------------------------------------------
// A tiling class.
//---------------------------------------------------------------------
public class Tiling {

    public int a, b, c, N, S, T;   // the usual parameters of the hexagon
    public int[][] particles; // an N x T array that contains the position of
    // all particles

    //---------------------------------------------------------------------
    // Create standard tiling (all paths go right and then up) of a hexagon
    // of side lengths a, b and c
    //---------------------------------------------------------------------
    public Tiling(int a, int b, int c) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.N = a;
        this.S = c;
        this.T = b + c;
        this.particles = new int[N][T + 1];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= T; j++) {
                if (j < b) {
                    particles[i][j] = i;
                } else {
                    particles[i][j] = i + (j - b);
                }

            }
        }

    }

    //---------------------------------------------------------------------
    // Return the particle locations
    //---------------------------------------------------------------------
    public int[][] getParticles() {
        return particles;
    }
    
    //---------------------------------------------------------------------
    // Validate the tiling
    //---------------------------------------------------------------------
    public boolean verify() {
        boolean check = true;
        for (int t = 1; t < (T + 1); t++) {
            for (int i = 0; i < N; i++) {
                if (particles[i][t] < particles[i][t - 1]) {
                    check = false;
                }
                if (particles[i][t] > (particles[i][t - 1] + 2)) {
                    check = false;
                }
            }
        }

        for (int i = 0; i < N; i++) {
            if (particles[i][0] != i) {
                check = false;
            }
            if (particles[i][T] != (i + S)) {
                System.out.println(particles[i][T] + " " + (i + S));
                check = false;
            }
        }
        return check;
    }


    //---------------------------------------------------------------------
    // Draw the tiling with squares and 45 degree sloped lines
    //---------------------------------------------------------------------
    public void draw45() {

        int scale = Math.max(T, N + S + 2);
        StdDraw.setXscale(-1, scale);
        StdDraw.setYscale(-1, scale);

        double b = 0.001;

        StdDraw.setPenRadius(b);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < T; j++) {
                double[] x = new double[4];
                double[] y = new double[4];
                int p = particles[i][j];
                // draw a "square" lozenge
                if (p == particles[i][j + 1]) {
                    x[0] = j;
                    y[0] = p;
                    x[1] = j + 1;
                    y[1] = p;
                    x[2] = j + 1;
                    y[2] = p + 1;
                    x[3] = j;
                    y[3] = p + 1;
                }
                // draw an "up" lozenge
                else {
                    x[0] = j;
                    y[0] = p;
                    x[1] = j;
                    y[1] = p + 1;
                    x[2] = j + 1;
                    y[2] = p + 2;
                    x[3] = j + 1;
                    y[3] = p + 1;
                }
                StdDraw.polygon(x, y);
            }
        }

        for (int j = 1; j < T - S; j++) {
            for (int k = 0; k < Math.min(N + j, S + N); k++) {
                boolean check = true;
                for (int i = 0; i < N; i++) {
                    if (particles[i][j] == k) {
                        check = false;
                    }
                }

                // draw a "flat" lozeange
                if (check) {
                    double[] x = new double[4];
                    double[] y = new double[4];
                    x[0] = j;
                    y[0] = k;
                    x[1] = j + 1;
                    y[1] = k + 1;
                    x[2] = j;
                    y[2] = k + 1;
                    x[3] = j - 1;
                    y[3] = k;
                    StdDraw.polygon(x, y);
                }
            }
        }

        for (int j = T - S; j < T; j++) {
            for (int k = j + S - T; k < Math.min(N + j, S + N); k++) {
                boolean check = true;
                for (int i = 0; i < N; i++) {
                    if (particles[i][j] == k) {
                        check = false;
                    }
                }
                // draw a "flat" lozeange
                if (check) {
                    double[] x = new double[4];
                    double[] y = new double[4];
                    x[0] = j;
                    y[0] = k;
                    x[1] = j + 1;
                    y[1] = k + 1;
                    x[2] = j;
                    y[2] = k + 1;
                    x[3] = j - 1;
                    y[3] = k;
                    StdDraw.polygon(x, y);
                }
            }
        }

    }

    //---------------------------------------------------------------------
    // Draw the tiling with squares and 45 degree sloped lines
    //---------------------------------------------------------------------
    public void drawLines() {

        StdDraw.setXscale(-1, T + 2);
        StdDraw.setYscale(-1, T + 2);
        double a = 0.01;
        double b = 0.001;

        StdDraw.setPenRadius(b);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < T; j++) {
                double[] x = new double[4];
                double[] y = new double[4];
                int p = particles[i][j];
                // draw a "square" lozenge
                if (p == particles[i][j + 1]) {
                    x[0] = j;
                    y[0] = p;
                    x[1] = j + 1;
                    y[1] = p;
                    x[2] = j + 1;
                    y[2] = p + 1;
                    x[3] = j;
                    y[3] = p + 1;

                    StdDraw.setPenRadius(a);
                    StdDraw.line(x[0], y[0] + 0.5, x[1], y[0] + 0.5);
                    StdDraw.setPenRadius(b);

                }
                // draw an "up" lozenge
                else {


                    x[0] = j;
                    y[0] = p;
                    x[1] = j;
                    y[1] = p + 1;
                    x[2] = j + 1;
                    y[2] = p + 2;
                    x[3] = j + 1;
                    y[3] = p + 1;

                    StdDraw.setPenRadius(a);
                    StdDraw.line(x[0], y[0] + 0.5, x[2], y[1] + 0.5);
                    StdDraw.setPenRadius(b);

                }
                StdDraw.polygon(x, y);
            }
        }

        for (int j = 1; j < T - S; j++) {
            for (int k = 0; k < Math.min(N + j, S + N); k++) {
                boolean check = true;
                for (int i = 0; i < N; i++) {
                    if (particles[i][j] == k) {
                        check = false;
                    }
                }

                // draw a "flat" lozenge
                if (check) {
                    double[] x = new double[4];
                    double[] y = new double[4];
                    x[0] = j;
                    y[0] = k;
                    x[1] = j + 1;
                    y[1] = k + 1;
                    x[2] = j;
                    y[2] = k + 1;
                    x[3] = j - 1;
                    y[3] = k;
                    StdDraw.polygon(x, y);
                }
            }
        }

        for (int j = T - S; j < T; j++) {
            for (int k = j + S - T; k < Math.min(N + j, S + N); k++) {
                boolean check = true;
                for (int i = 0; i < N; i++) {
                    if (particles[i][j] == k) {
                        check = false;
                    }
                }
                // draw a "flat" lozenge
                if (check) {
                    double[] x = new double[4];
                    double[] y = new double[4];
                    x[0] = j;
                    y[0] = k;
                    x[1] = j + 1;
                    y[1] = k + 1;
                    x[2] = j;
                    y[2] = k + 1;
                    x[3] = j - 1;
                    y[3] = k;
                    StdDraw.polygon(x, y);
                }
            }
        }

    }

    //---------------------------------------------------------------------
    // Return the min and max coordinates for holes on each vertical level
    //---------------------------------------------------------------------
    public int[][] boundary() {
        int[][] coord = new int[T][2];

        for (int j = 1; j < T; j++) {
            int min = Math.min(N + j, N + S);
            int max = 0;
            boolean check = false;
            for (int i = 0; i < (N - 1); i++) {
                if (particles[i + 1][j] != (particles[i][j] + 1)) {
                    if (min > particles[i][j]) {
                        min = particles[i][j];
                        check = true;
                    }
                    if (max < particles[i + 1][j]) {
                        max = particles[i + 1][j];
                    }
                }

            }
            if (!check) {
                StdOut.println(j);
            }
            coord[j][0] = min;
            coord[j][1] = max;
        }
        return coord;
    }


    //---------------------------------------------------------------------
    // Print the tiling
    //---------------------------------------------------------------------
    public void print() {
        StdOut.println(a + " " + b + " " + c);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < (T + 1); j++) {
                StdOut.print(particles[N - i - 1][j] + " ");
            }
            StdOut.println();
        }
    }

    //---------------------------------------------------------------------
    // Converts the tiling to a string representation
    //---------------------------------------------------------------------
    public static void main(String[] args) {
        int a = 3;
        int b = 3;
        int c = 3;
        Tiling tile = new Tiling(a, b, c);
        int[][] p = tile.particles;
        p[0][0] = 0;
        p[0][1] = 0;
        p[0][2] = 0;
        p[0][3] = 1;
        p[0][4] = 1;
        p[0][5] = 2;
        p[0][6] = 3;
        p[1][0] = 1;
        p[1][1] = 1;
        p[1][2] = 2;
        p[1][3] = 3;
        p[1][4] = 4;
        p[1][5] = 4;
        p[1][6] = 4;
        p[2][0] = 2;
        p[2][1] = 3;
        p[2][2] = 3;
        p[2][3] = 4;
        p[2][4] = 5;
        p[2][5] = 5;
        p[2][6] = 5;

        tile.draw45();
        tile.print();
    }
}
