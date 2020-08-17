import java.lang.Math;
public class Main {
    
    public static void main(String args[]) {
        FFT fft=new FFT(8);
        double[] a=new double[8];
        a[0]=1.0;
        System.out.println( fft.power_spector(a) );
    }
}