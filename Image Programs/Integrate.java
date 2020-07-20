//---------------------------------------------------------------------
// Integrate class: Contains two functions for integrating covariance
// on pairs of vertical slices for the qHahn case.
//---------------------------------------------------------------------

public class Integrate{
    
    //---------------------------------------------------------------------
    // Sample functions to be integrated: polynomials and cosine
    //---------------------------------------------------------------------
    public double fun(double x, int n){
        double f = 1;
        for (int i = 0; i < n; i++){
            f = f*x;
        }
        return f;
    }
    
    public double cos(double x, int n){
        
        return n*Math.cos(n*x);
    }
    
    
    //---------------------------------------------------------------------
    // Integrator based on Riemann sums
    //---------------------------------------------------------------------
    public void Int2(){
        Integrate I = new Integrate();
        double q0 = 0.998;
        double q = I.fun(q0, 500);
        double S = I.fun(1/q0, 500);
        double T = I.fun(1/q0, 1000);
        double k = 0.1;
        
        double x1 = I.fun(1/q0, 500);
        double x2 = I.fun(1/q0, 250);
        
        double dC1 = q*q*16*(1-q)*(1-k*q)*(x1 - 1)*(T - q)*(1-k*T*q)*(T - x1)*(S - 1)*(T-S);
        double aC1 = q*q*(T - 1)*(T - 1);
        double bC1 = 2*(T-1)*(S*x1 - T +(2 - S - x1)*q);
        bC1 = bC1 - 4*(T - q)*(S - 1)*(x1 - 1);
        bC1 = bC1*q - 2*q*k*(S*T*x1*q*q + S*T*T*q - 2*S*T*x1*q -2*S*T*q*q);
        bC1 = bC1 - 2*q*k*(S*x1*q*q + T*T*x1*q + T*T*q*q - 2*T*x1*q*q);
        bC1 = bC1 - 2*q*k*(S*T*q - 2*T*T*q + T*x1*q + T*q*q);
    
        double dC2 = q*q*16*(1-q)*(1-k*q)*(x2 - 1)*(T - q)*(1-k*T*q)*(T - x2)*(S - 1)*(T-S);
        double aC2 = q*q*(T - 1)*(T - 1);
        double bC2 = 2*(T-1)*(S*x2 - T +(2 - S - x2)*q);
        bC2 = bC2 - 4*(T - q)*(S - 1)*(x2 - 1);
        bC2 = bC2*q - 2*q*k*(S*T*x2*q*q + S*T*T*q - 2*S*T*x2*q -2*S*T*q*q);
        bC2 = bC2 - 2*q*k*(S*x2*q*q + T*T*x2*q + T*T*q*q - 2*T*x2*q*q);
        bC2 = bC2 - 2*q*k*(S*T*q - 2*T*T*q + T*x2*q + T*q*q);
        
        double a1 = ((-1)*bC1 - Math.sqrt(dC1))/(2*aC1);
        double b1 = ((-1)*bC1 + Math.sqrt(dC1))/(2*aC1);
        double a2 = ((-1)*bC2 - Math.sqrt(dC2))/(2*aC2);
        double b2 = ((-1)*bC2 + Math.sqrt(dC2))/(2*aC2);
        
        StdOut.println(a1 + " " + b1);
        int deg = 5;
        int iter = 1000;
        double h1 = (b1 - a1) / iter;     // step sizes
        double h2 = (b2 - a2) / iter;     // step sizes

        for (int i1 = 0; i1 < deg; i1++){
            for (int i2 = 0; i2 < deg; i2++){
                double sum = 0;
                for (int i = 1; i < (iter-1); i++){
                    for (int j = 1; j < (iter-1); j++){
                        double y1 = a1 + h1*i;              // these are the coordinates Y + S*k*X/Y
                        double y2 = a2 + h2*j;
                        y1 = (y1 + Math.sqrt(y1*y1 - 4*S*k*x1))/2; //  these are the true Y coordinates
                        y2 = (y2 + Math.sqrt(y2*y2 - 4*S*k*x2))/2;
                        

                        double aM1 = q*(y1-1)*(1 - S*k*x1/y1);
                        double aM2 = q*(y2-1)*(1 - S*k*x2/y2);
                        
                        double bM1 = q*(T-1)*y1 + S*x1 - S*q - x1*q - T + 2*q;
                        bM1 = bM1 -q*(2*S*T*x1 - S*x1*q - S*T - T*x1 + T*q)*k - q*k*(S*x1 - S*T*x1)/y1;
                        double bM2 = q*(T-1)*y2 + S*x2 - S*q - x2*q - T + 2*q;
                        bM2 = bM2 -q*(2*S*T*x2 - S*x2*q - S*T - T*x2 + T*q)*k - q*k*(S*x2 - S*T*x2)/y2;
                        
                        double cM1 = (T-q)*(S-1)*(x1-1)*(1-T*k*q);
                        double cM2 = (T-q)*(S-1)*(x2-1)*(1-T*k*q);
                        
                        double dM1 = Math.sqrt(4*aM1*cM1 - bM1*bM1);
                        double dM2 = Math.sqrt(4*aM2*cM2 - bM2*bM2);

                        double rNum = -bM1/(2*aM1) + bM2/(2*aM2);
                        double rDen = rNum;
                        double iNum = dM1/(2*aM1) - dM2/(2*aM2);
                        double iDen = dM1/(2*aM1) + dM2/(2*aM2);
                        
                        double fun1 = I.fun(y1, i1 );
                        double fun2 = I.fun(y2, i2 );
                        
                        sum = sum + 0.5*(Math.log(rNum*rNum + iNum*iNum) - Math.log(rDen*rDen + iDen*iDen))*fun1*fun2;
                    }

                }
                
                sum = -sum*h1*h2/(2*Math.PI*Math.PI);
                System.out.print(sum + " ");
                
            }
            System.out.println();
        }
    }
   
    
    //---------------------------------------------------------------------
    // Integrator based on random sampling
    //---------------------------------------------------------------------
    public void Int(){
        Integrate I = new Integrate();
        
        double q0 = 0.998;
        double q = I.fun(q0, 500);
        double S = I.fun(1/q0, 500);
        double T = I.fun(1/q0, 1000);
        
        double x1 = I.fun(1/q0, 250);
        double x2 = I.fun(1/q0, 500);
        
        double aC1 = q*q*(T-1)*(T-1);
        double aC2 = aC1;
        
        double bC1 = 2*((x1 -1)*(S-q) - q*(T-1) + (S-T)*(1-q))*q*(T-1) - 4*q*(T-q)*(S-1)*(x1-1);
        double bC2 = 2*((x2 -1)*(S-q) - q*(T-1) + (S-T)*(1-q))*q*(T-1) - 4*q*(T-q)*(S-1)*(x2-1);
        
        double cC1 = (x1 - 1)*(S-q) - q*(T-1) + (S-T)*(1-q);
        cC1 = cC1*cC1 + 4*q*(T-q)*(S-1)*(x1-1);
        double cC2 = (x2 - 1)*(S-q) - q*(T-1) + (S-T)*(1-q);
        cC2 = cC2*cC2 + 4*q*(T-q)*(S-1)*(x2-1);
        
        double dC1 = bC1*bC1 - 4*aC1*cC1;
        double dC2 =  bC2*bC2 - 4*aC2*cC2;
        
        double a1 = ((-1)*bC1 - Math.sqrt(dC1))/(2*aC1);
        double b1 = ((-1)*bC1 + Math.sqrt(dC1))/(2*aC1);
        double a2 = ((-1)*bC2 - Math.sqrt(dC2))/(2*aC2);
        double b2 = ((-1)*bC2 + Math.sqrt(dC2))/(2*aC2);
        
        StdOut.println(a1 + " " + b1);
        int iter = 1000000;
        
        int deg = 5;
        
        for (int i1 = 0; i1 < deg; i1++){
            for (int i2 = 0; i2 < deg; i2++){
                double sum = 0;
                for (int i = 0; i < iter; i++){
                    double y1 = a1 + (b1 - a1)*Math.random();
                    double y2 = a2 + (b2 - a2)*Math.random();
                    
                    double aM1 = q*(y1-1);
                    double aM2 = q*(y2-1);
                    
                    double bM1 = (x1 - 1)*(S-q) + q*(T-1)*(y1 - 1) + (S-T)*(1-q);
                    double bM2 = (x2 - 1)*(S-q) + q*(T-1)*(y2 - 1) + (S-T)*(1-q);
                    
                    double cM1 = (T-q)*(S-1)*(x1-1);
                    double cM2 = (T-q)*(S-1)*(x2-1);
                    
                    double dM1 = Math.sqrt(4*aM1*cM1 - bM1*bM1);
                    double dM2 = Math.sqrt(4*aM2*cM2 - bM2*bM2);
                    
                    double rNum = -bM1/(2*aM1) + bM2/(2*aM2);
                    double rDen = rNum;
                    double iNum = dM1/(2*aM1) - dM2/(2*aM2);
                    double iDen = dM1/(2*aM1) + dM2/(2*aM2);
                    
                    double fun1 = I.fun(y1, i1 );
                    double fun2 = I.fun(y2, i2 );
                    
                    double f1 = Math.sqrt((y1-a1)*(b1-y2));
                    double f2 = Math.sqrt((y2-a1)*(b1-y1));
                       

                    sum = sum + 0.5*(Math.log(rNum*rNum + iNum*iNum) - Math.log(rDen*rDen + iDen*iDen))*fun1*fun2;
                    
                }
                
                sum = -sum*(b1-a1)*(b2-a2)/(2*Math.PI*Math.PI*iter);
                System.out.print(sum + " ");
            }
            System.out.println();
        }
    }
    
    
    
    public static void main (String[] args){
        Integrate I = new Integrate();
        I.Int2();
    }
    
    
}