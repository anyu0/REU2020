import java.util.Arrays;

//---------------------------------------------------------------------
// testFluctuations: Contains functions that can be used to study
// the limit shape (LS below) the law of large numbers (LLN below)
// and the central limit theorem (CLT 1-4 below) for random lozenge 
// tilings with the q-Racah distribution
//---------------------------------------------------------------------
public class testFluctuations{
    
    
    // constructor
    public testFluctuations(){
    }
    
    
    //---------------------------------------------------------------------
    // Some possible test functions
    //---------------------------------------------------------------------
    public double poly(double x, int n){
        double f = 1;
        for (int i = 0; i < n; i++){
            f = f*x;
        }
        return f/n;
    }
    
    public double sine(double x, int n){
        return Math.sin(n*x);
    }
    
    //---------------------------------------------------------------------
    // Sample limit shape vs actual limit shape
    //---------------------------------------------------------------------
    public void LS(){
        double q = StdIn.readDouble();
        double k = StdIn.readDouble();
        int N = StdIn.readInt();
        int T = StdIn.readInt();
        int S = StdIn.readInt();
        int iter = 100;
        int M1 = 5000;
        
        double[] qM = new double[4*M1 + 1]; 
        qM[2*M1] = 1;
        for (int i = 1; i < 2*M1; i++){          
            qM[i + 2*M1] = qM[i-1 + 2*M1]*q;
            qM[-i + 2*M1] = qM[-i + 1 + 2*M1]/q;        
        }
        
        int M = Math.max(N+S,T);
        
        StdDraw.setXscale(1, qM[2*M1-M]*(1+k));
        StdDraw.setYscale(1, qM[2*M1-M]*(1+k));
        
        // sample iter many tilings and the top and bottom hole position
        double[][][] pts = new double[iter][T+1][4];
        for (int i = 0; i < iter; i++){
            Tiling tile = new Tiling(N, T-S,S);
            for (int i1 = 0; i1 < N; i1++){
                for (int j1 = 0; j1 < (T+1); j1++){
                    tile.particles[i1][j1] = StdIn.readInt();
                }
            }
            int[][] coord = tile.boundary();
            for (int j = 1; j < T; j++){
                double x = qM[coord[j][0] + 2*M1];
                double y = qM[coord[j][1] + 2*M1];
                
                x = 1/x + x*k*qM[2*M1 - S - j];
                y = 1/y + y*k*qM[2*M1 - S - j];
                                
                pts[i][j][0] = x;
                pts[i][j][1] = y;
            }
        }
        
        double[] avx = new double[T+1];
        double[] avy = new double[T+1]; 
        for (int i = 0; i < iter; i++){
            for (int j = 1; j < T; j++){
                avx[j] = avx[j] + pts[i][j][0];
                avy[j] = avy[j] + pts[i][j][1];
            }
        }
        
        // draw empirical limit shape
        for (int j = 1; j < T; j++){
            avx[j] = avx[j]/iter;
            avy[j] = avy[j]/iter;
            StdDraw.filledCircle(qM[2*M1 -j ], avx[j] , qM[2*M1 -M]*0.003);
            StdDraw.filledCircle(qM[2*M1-j], avy[j] , qM[2*M1 -M]*0.003);
        }
        
        // draw theoretical limit shape
        StdDraw.setPenColor(StdDraw.BOOK_RED);
        double T1 = qM[2*M1-T];
        double S1 = qM[2*M1-S];
        double q1 = qM[2*M1 + N];
        
        for (int j = 1; j < T; j++){
            double x1 = qM[2*M1-j];
            double D = q1*Math.sqrt(16*(1-q1)*(1-k*q1)*(x1 - 1)*(T1 - q1)*(1-k*T1*q1)*(T1 - x1)*(S1 - 1)*(T1-S1));
            double a1 = q1*q1*(T1 - 1)*(T1 - 1);
            double b1 = 2*(T1-1)*(S1*x1 - T1 +(2 - S1 - x1)*q1);
            b1 = b1 - 4*(T1 - q1)*(S1 - 1)*(x1 - 1);
            b1 = b1*q1 - 2*q1*k*(S1*T1*x1*q1*q1 + S1*T1*T1*q1 - 2*S1*T1*x1*q1 -2*S1*T1*q1*q1);
            b1 = b1 - 2*q1*k*(S1*x1*q1*q1 + T1*T1*x1*q1 + T1*T1*q1*q1 - 2*T1*x1*q1*q1);
            b1 = b1 - 2*q1*k*(S1*T1*q1 - 2*T1*T1*q1 + T1*x1*q1 + T1*q1*q1);

            double top = (-b1 + D)/(2*a1);
            double bot = (-b1 - D)/(2*a1);
            
            StdDraw.filledCircle(qM[2*M1-j], top ,qM[2*M1-M]*0.003);
            StdDraw.filledCircle(qM[2*M1-j], bot , qM[2*M1-M]*0.003);
        }
        
        
    }
    
    //---------------------------------------------------------------------
    // Calculate the covariance for a given pair of slices
    //--------------------------------------------------------------------- 
    public void CLT(){
        double q = StdIn.readDouble();
        double k = StdIn.readDouble();
        int N = StdIn.readInt();
        int T = StdIn.readInt();
        int S = StdIn.readInt();
        int iter = 500;
        int M = 5000;
        int deg = 5;
        
        double[][] fun1 = new double[iter][deg];
        double[][] fun2 = new double[iter][deg];
        double[] qM = new double[M];
        
        qM[0] = 1;
        for (int i = 1; i < M; i++){
            qM[i] = qM[i-1]/q;
        }
        
        
        // Read the particle configurations
        for (int i = 0; i < iter; i++){
            int[][] diag = new int[N][T+1];
            for (int j = 0; j < N; j++){
                for (int t = 0; t < (T+1); t++){
                    diag[j][t] = StdIn.readInt();
                }
            }
            
            testFluctuations f = new testFluctuations();
            
            // the two vertical slices
            int t1 = T/2;
            int t2 = T/4;
            // evaluate the functions on the given two slices
            for (int u = 0; u < deg; u++){
                double sum1 = 0;
                double sum2 = 0;
                for (int j = 0; j < N; j++){ 
                    sum1 = sum1 + f.poly(qM[diag[j][t1]] + k*qM[j+S]/qM[diag[j][t1]], u + 1);
                    sum2 = sum2 + f.poly(qM[diag[j][t2]]+ k*qM[j+S]/qM[diag[j][t2]], u + 1);
                }
                fun1[i][u] = sum1/N;
                fun2[i][u] = sum2/N;
            }
            
        }
        
        // calculate and print the covariances
        for (int u = 0; u < deg; u++){
            for (int v = 0; v < deg; v++){
                double cov = 0;
                double sum1 = 0;
                double sum2 = 0;
                double prod = 0;
                for (int i = 0; i < iter; i++){
                    sum1 = sum1 + fun1[i][u];
                    sum2 = sum2 + fun2[i][v];
                }
                
                sum1 = sum1/iter;
                sum2 = sum2/iter;
                for (int i = 0; i < iter; i++){
                    prod = prod + (fun1[i][u] - sum1)*(fun2[i][v] - sum2);
                }
                cov = N*N*prod/iter;
                StdOut.print("Case:" + u + " " + v + " Cov: " + cov + " ");
            }
            StdOut.println();
        }

    }
    
    //---------------------------------------------------------------------
    // A function that tests the law of large numbers - it takes samples,
    // calculates the empirical density and compares it with the theore-
    // tical one.
    //--------------------------------------------------------------------- 
    
    public void LLN(){
        double q = StdIn.readDouble();
        double k = StdIn.readDouble();
        int N = StdIn.readInt();
        int T = StdIn.readInt();
        int S = StdIn.readInt();
        int iter = 500;
        int M = 5000;
        
        double[] qM = new double[M]; 
        qM[0] = 1;
        for (int i = 1; i < M; i++){
            qM[i] = qM[i-1]/q;
        }
        
        // the vertical slice on which we calculate the density
        int t1 = T/2;
        
        // store the frequencies of each particle position
        double[] freq = new double[N + S];
        for (int i = 0; i < iter; i++){
            int[][] diag = new int[N][T+1];
            for (int j = 0; j < N; j++){
                for (int t = 0; t < (T+1); t++){
                    diag[j][t] = StdIn.readInt();
                }
            }
            
            testFluctuations f = new testFluctuations();
            

            for (int j = 0; j < N; j++){
                freq[diag[j][t1]]++;
            }                  
        }
        
        // normalize the frequency and put it in a histogram
        int r = 50;
        int bins = 2*r;
        double a = 1;
        double b = qM[S+N];
        double h = 1.0*(b-a)/bins;
        
        StdDraw.setXscale(a, b);
        StdDraw.setYscale(0,1);
        
        double[] eden = new double[bins];
        for (int i = 0; i < (N+S); i++){    
            for (int j = 0; j < bins; j++){
                double x = qM[i];
                if ((x <  a + (j + 1)*h) && ( x > a+ j*h)){
                    eden[j] = eden[j] + freq[i]/(h*N*iter);
                }
            }
        }
        for (int j = 0; j < bins; j++){
       
            double x1 =  a + j*h + h/2;
            double y1 = eden[j]/2;
            
            StdDraw.rectangle(x1, y1, h/2, y1); 
            StdOut.print(2*y1 + " ");
            
            // evaluate the theoretical density
            double T1 = qM[T];
            double S1 = qM[S];
            double N1 = qM[N];
            double x = x1;
            double t = qM[t1];  
            double A = (1 - S1*N1/x)*(1-t*N1/x)*(1-k*T1/x)*(1-k*t*S1/x);
            double B = N1*N1*T1*(1 - 1/x)*(1-t*S1/(x*T1))*(1-k*t/(N1*x))*(1 - k*S1/(N1*x));
            double phi = (N1 - 1)*(1 - T1*N1)*(1-k*t*S1/(x*x))*(1-k*t*S1/(x*x)) + A+ B;
            double den = 0;
            if ((A*B) == 0){
                if (phi > 0) den = 0;
                if (phi < 0) den = Math.PI;
            }
            else if (A*B > 0){
                phi = phi/(2*Math.sqrt(A*B));
                if (phi > 1) den = 0;
                else if (phi < -1) den = Math.PI;
                else den = Math.acos(phi);
            }
            else{
                den = 0;
            }
            
            // den is the theoretical density
            den = den/(Math.PI*x);
            StdOut.println(den + " ");
            
            
            StdDraw.setPenColor(StdDraw.BOOK_RED);
            double y2 = den;
            StdDraw.filledEllipse(x1, y2,(b-a)*0.01, 0.01);
            StdDraw.setPenColor();
        }
        
        StdOut.println();
 
    }
    
//*************************************************************************************
// Some additional functions that I've written for easing up calculations
// they are not important.
//*************************************************************************************
    
    public void CLT2(){
        double q = StdIn.readDouble();
        double k = StdIn.readDouble();
        int N = StdIn.readInt();
        int T = StdIn.readInt();
        int S = StdIn.readInt();
        int iter = 5000;
        int M = 5000;
        int deg = 5;
        
        double[][] fun1 = new double[iter][deg];
        double[][] fun2 = new double[iter][deg];
        double[] qM = new double[M];
        
        qM[0] = 1;
        for (int i = 1; i < M; i++){
            qM[i] = qM[i-1]/q;
        }
        
        
        for (int i = 0; i < iter; i++){
            int[][] diag = new int[N][T+1];
            for (int j = 0; j < N; j++){
                for (int t = 0; t < (T+1); t++){
                    diag[j][t] = StdIn.readInt();
                }
            }
            
            testFluctuations f = new testFluctuations();
            
            for (int t1 = 0; t1 < (T+1); t1 = t1 + 10){
                for (int u = 0; u < deg; u++){
                    double sum1 = 0;
                    for (int j = 0; j < N; j++){
                        sum1 = sum1 + f.poly(qM[diag[j][t1]], u + 1);
                    }
                    StdOut.print(sum1/N + " ");
                }
                StdOut.println();
            }
        }        
        
    }
   
    public void CLT3(){
        int iter = 500;
        int deg = 5;
        int t1 = 40;
        int t2 = 60;
        int total = 101;
        int N = 500;
        
        double[][][] data = new double[iter][deg][total];
        
        for (int i = 0; i < iter; i++){
            for (int t = 0; t < total; t++){
                for (int u = 0; u < deg; u++){
                    data[i][u][t] = StdIn.readDouble();
                }
            }
        }    
        
        deg = 1;
        for (int u = 0; u < deg; u++){
            for (int v = 0; v < deg; v++){
                double cov = 0;
                double sum1 = 0;
                double sum2 = 0;
                double prod = 0;
                for (int i = 0; i < iter; i++){
                    sum1 = sum1 + data[i][u][t1];
                    sum2 = sum2 + data[i][v][t2];
                    //prod = prod + data[i][u][t1]*data[i][v][t2];
                }
                
                sum1 = sum1/iter;
                sum2 = sum2/iter;
                for (int i = 0; i < iter; i++){
                    StdOut.println(iter*(data[i][u][t1] - sum1)*(data[i][u][t1] - sum2));
                }
                
            }
            StdOut.println();
        }
        
    }
    
    public void CLT4(){
        int iter = 5000;
        int t1 = 4;
        int total = 41;
        int N = 200;
        
        int deg = 5;
        int pow = 0;
        double[][][] data = new double[iter][deg][total];
        
        for (int i = 0; i < iter; i++){
            for (int t = 0; t < total; t++){
                for (int u = 0; u < deg; u++){
                    data[i][u][t] = StdIn.readDouble();
                }
            }
        }    
        
        
        for (int u = 0; u < deg; u++){
            double var = 0;
            double mean = 0;
            for (int i = 0; i < iter; i++){
                mean = mean + data[i][u][t1];
            }
            
            mean = mean/iter;
            for (int i = 0; i < iter; i++){
                var = var + (data[i][u][t1] - mean)*(data[i][u][t1] - mean);
            }
            
            var = var/iter;
            
            StdOut.println("Variance = " + N*N*var);
        }
        
        double var = 0;
        double mean = 0;
        for (int i = 0; i < iter; i++){
            mean = mean + data[i][pow][t1];
        }
        
        mean = mean/iter;
        for (int i = 0; i < iter; i++){
            var = var + (data[i][pow][t1] - mean)*(data[i][pow][t1] - mean);
        }
        var = var/iter;
        double sdev = Math.sqrt(var);
        double[] result = new double[iter];
        for (int i = 0; i < iter; i++){
            result[i] = (data[i][pow][t1] - mean)/sdev;
        }   
        
        Arrays.sort(result);
        for (int i = 0; i < iter; i++){
            StdOut.println(result[i]);
        }   
        
        int r = 20;
        int bins = 2*r;
        double a = 3;
        double h = 1.0*2*a/bins;
        
        StdDraw.setXscale(-a, a);
        StdDraw.setYscale(0,1);
        
        double[] eden = new double[bins];
        for (int i = 0; i < iter; i++){    
            for (int j = 0; j < bins; j++){
                double x = result[i];
                if ((x <  (j-r + 1)*h) && ( x >  (j - r)*h)){
                    eden[j] = eden[j] + 1.0/(h*iter);
                }
            }
        }
        for (int j = 0; j < bins; j++){
            
            double x1 =  (j - r)*h + h/2;
            double y1 = eden[j]/2;
                        
            StdDraw.rectangle(x1, y1, h/2, y1); 
            double y2 = Math.exp(-x1*x1/2)/Math.sqrt(2*Math.PI);
            StdDraw.filledCircle(x1, y2, 0.01);
        }
        
        
    }
    
    //--------------------------------------------------------------------------
    // The main method: Do what you want here
    //--------------------------------------------------------------------------
    public static void main(String[] args){
        testFluctuations f = new testFluctuations();
        f.CLT();  
    }
}