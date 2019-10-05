package cs3010;
import java.util.*;
import java.io.*;

public class polRoot {

	public static int iterCounter;
	public static int defaultIter = 10000;
	public static final float machineEPS = 2E-23F;
	
	public static void main(String[] args) {

		float poly[] = {1, 2, 10, -20};
		System.out.println(Bisection(poly, 0, 2, defaultIter, machineEPS));
		System.out.println(Newton(poly, DeriveFunction(poly), 2,  defaultIter, machineEPS, 0));
		System.out.println(Secant(poly, 2,4, defaultIter, machineEPS)); 
		System.exit(0);

	}
	
	//Solves polynomials via Horner's algorithm
	static float Function(float[] coeff, float x) {
		float result = coeff[0];
		
		for(int i = 1; i < coeff.length ; i++) {
			result = x * result + coeff[i];
		}
		
		return result;
	}
	static float[] DeriveFunction(float[] coeff) {
		float derivative[] = new float[coeff.length-1];
		
		for(int i = 0; i < derivative.length; i++) {
			derivative[i] = coeff[i] * (coeff.length - 1 - i);
		}
		
		//Decrease array size if coefficient reults in 0 from deriving
		if(derivative[derivative.length-1] == 0) {
			
			derivative = Arrays.copyOf(derivative, derivative.length-1);
		}
		
		return derivative;
	}
	static float Bisection(float[] f, float a, float b, int maxIter, float eps) {

		float fa = Function(f, a);
		float fb = Function(f, b);
		
		if( fa * fb >= 0) {
			System.out.println("Inadequate values for a and b.");
			return -1;
		}
		
		float error = b - a;
		float c = 0;
		for(int i = 0; i < maxIter; i++) {
			error = error / 2;
			c = a + error;
			float fc = Function(f, c);
			
			iterCounter++;
			if( Math.abs(error) < eps || fc == 0) {
				System.out.println("Algorithm has converged after " + i + " iterations!");
				return c;
			} //end if
			
			if( fa * fc < 0) {
				b = c;
				fb = fc;
			}
			else {
				a = c;
				fa = fc;
			} //end if
		
		} //end for
		
		System.out.println("Max iterations reached without convergence...");
		return c;
	} //end Bisection method
	static float Newton(float[] f, float[] df, float x, int maxIter, float eps, float delta) {
		float fx = Function(f, x);
		for(int i = 0; i < maxIter; i++) {
			float fd = Function(df, x);
			if ( Math.abs(fd) < delta) {
				System.out.println("Small slope!");
				return x;
			} //end if
			
			float d = (float)(fx / fd);
			x = x - d;
			fx = Function(f, x);
			
			if(Math.abs(d) < eps) {
				System.out.println("Algorithm has converged after " + i + " iterations!");
				return x;
			} //end if
		} //end for
		
		System.out.println("Max iterations reached without convergence...");
		return x;
	} //end Newton Method
	static float Secant(float[] f, float a, float b, int maxIter, float eps) {
		float fa = Function(f, a);
		float fb = Function(f, b);
		
		if(Math.abs(fa) > Math.abs(fb)) {
			// swap(a, b)
			float temp = a;
			a = b;
			b = temp;
			// sap (fa, fb)
			float temp1 = fa;
			fa = fb;
			fb = temp1;
		} //end if
		
		for(int i = 0; i < maxIter; i++) {
			if(Math.abs(fa) > Math.abs(fb)) {
				// swap(a, b)
				float temp = a;
				a = b;
				b = temp;
				// sap (fa, fb)
				float temp1 = fa;
				fa = fb;
				fb = temp1;
			} //end if
			
			float d = (b - a) / (fb - fa);
			b = a;
			fb = fa;
			d = d * fa;
			
			if(Math.abs(d) < eps) {
				System.out.println("Algorithm has converged after " + i + " iterations!");
				return a;
			} //end if
			
			a = a - d;
			fa = Function(f, a);
		} //end for
		
		System.out.print("Maximum number of itertions reached!");
		return a;
	} //end Secant method
	static float Hybrid(float[] f, float[] df, float a, float b, int maxIter, float eps, float delta) {
		
	}

}
