public class DrawAvoidingBridges {


    public static void main(String[] args) {

        double q = 0.999999;
        double k = 0;
        int N = 20;
        int T = 20;
        int S = 10;
        SimulateQR qr = new SimulateQR(q, k, N, T, S);
        Tiling tile = qr.sample();
        int[][] pos = tile.getParticles();

        double[][] bb2 = new double[N][T];
        double slope = 1.0 * S / T;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < T; j++) {
                bb2[i][j] = pos[i][j] - slope * j - i;
            }
        }
        tile.drawLines();
        // draw the trajectories
        /*StdDraw.setXscale(-100, 600);
        StdDraw.setYscale(-100, 100);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < T - 1; j++) {
                StdDraw.line(j, bb2[i][j], j + 1, bb2[i][j + 1]);
            }
        }*/
    }


}
