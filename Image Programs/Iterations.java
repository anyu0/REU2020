import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

//---------------------------------------------------------------------
// Iterations class: Sample q-Racah lozenge tilings efficiently
// using some paralel computing.
//---------------------------------------------------------------------
public class Iterations extends Thread{
    private int N, S, T, iter, index;
    private double q,k;

    
   // Constructor
    public Iterations(double q, double k, int N, int T, int S, int iter, int index) {
        this.N = N;
        this.S = S;
        this.T = T;
        this.q = q;
        this.k = k;
        this.iter = iter;
        this.index = index;
    }
    
   
    // the sampling method
    public void run() {
        File file = new File("result" + this.index + ".txt");
 
        try {
            SimulateQR qr = new SimulateQR(q,k,N,T,S);
            
            FileWriter fwt = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bwt = new BufferedWriter(fwt);
 
            for (int i = 0; i < iter; i++){
                Tiling tile = qr.sample();
                for (int j = 0; j < N; j++){
                    for (int t = 0; t < (T+1); t++){
                        bwt.write(tile.particles[j][t] + " ");
                    }
                    bwt.newLine();
                }
            }      
            
            bwt.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        
    }
    
    public static void main(String[] args) throws Exception{
        Stopwatch st = new Stopwatch();   // a watch to measure time efficiency
       
        // the parameters for the q-Racah distribution
        int N = 500;
        int S = 500;
        int T = 1000;
        double q = 0.998;
        double k = 0;
        
        int iter = 1;  // the number of iterations
        int nT = 1;    // number of parallel processes
                
        Iterations[] iT = new Iterations[nT];
        
        for (int i = 0; i < nT; i++){
            iT[i] = new Iterations(q,k,N,T,S,iter, i);
            iT[i].start();
        }
  
        
        try { 
            for (int i = 0; i < nT; i++){
                iT[i].join();
            }
        } catch(InterruptedException e){}

        //System.out.println(st.elapsedTime());
    }
    
}