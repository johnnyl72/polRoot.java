package cs3010;
import java.util.*;
import java.io.*;

public class polRoot {

	static int defaultIter = 10000;
	static final float machineEPS = 2E-23F;
	
	public static void main(String[] args) {

		float poly[] = {10, 2, -2};
		System.out.println(Bisection(poly, 0, 5, defaultIter, machineEPS));
		System.out.println(Newton(poly, DeriveFunction(poly), 3,  defaultIter, machineEPS, 0));
		//System.out.println(Secant(poly, 2,1, defaultIter, machineEPS)); 
		System.out.println(Hybrid(poly, 0, 5, defaultIter, machineEPS, 0));
		
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
		int iterCounter = 0;
		if( fa * fb >= 0) {
			System.out.println("Inadequate values for a and b.");
			return -1;
		}
		
		float error = b - a;
		float c = 0;
		for(int i = 0; i < maxIter; i++) {
			iterCounter++;
			error = error / 2;
			c = a + error;
			float fc = Function(f, c);
		
			
			if( Math.abs(error) < eps || fc == 0) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!");
				return c;
			} //end if
			
			
			if( fa * fc < 0) {
				b = c;
				fb = fc;
			}
			else {
				a = c;
				fa = fc;
			} //end if-else
		
		} //end for
		
		System.out.println("Max iterations reached without convergence...");
		return c;
	} //end Bisection method
	static float Newton(float[] f, float[] df, float x, int maxIter, float eps, float delta) {
		float fx = Function(f, x);
		int iterCounter = 0;
		for(int i = 0; i < maxIter; i++) {
			iterCounter++;
			float fd = Function(df, x);
			if ( Math.abs(fd) < delta) {
				System.out.println("Small slope!");
				return x;
			} //end if
			
			float d = (float)(fx / fd);
			x = x - d;
			fx = Function(f, x);
			
			if(Math.abs(d) < eps) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!");
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
		int iterCounter = 0;
		for(int i = 0; i < maxIter; i++) {
			iterCounter++;
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
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!");
				return a;
			} //end if
			
			a = a - d;
			fa = Function(f, a);
		} //end for
		
		System.out.print("Maximum number of itertions reached!");
		return a;
	} //end Secant method
	static float Hybrid(float[] f, float a, float b, int maxIter, float eps, float delta) {
		float fa = Function(f, a);
		float fb = Function(f, b);
		float c = 0;
		float fc = 0;
		int iterCounter = 0 ;
		
		for(int i = 0; i < maxIter; i++) {
		iterCounter++;
		if( fa * fb >= 0) {
			System.out.println("Inadequate values for a and b.");
			return -1;
		}
		
		float error = b - a;
		c = 0;
		error = error / 2;
		c = a + error;
		fc = Function(f, c);
		
		if( Math.abs(error) < eps || fc == 0) {
		
			System.out.println("Algorithm has converged after " + iterCounter + " iterations!!!");
			return c;
		} //end if
		else if(Math.abs(error) < 1) {
			break;
		} //end if
		if( fa * fc < 0) {
			b = c;
			fb = fc;
		}
		else {
			a = c;
			fa = fc;
		} //end if-else
		}
		//Newton's Part
		for(int j = 0; j < maxIter; j++) {
			iterCounter++;
			float fd = Function(DeriveFunction(f), c);
			if ( Math.abs(fd) < delta) {
				//Continue with bisection method instead by brekaing out Newton's for loop
				break; 
			} //end if
			
			float d = (float)(fc / fd);
			//New c
			c = c - d;
			fc = Function(f, c);
			
			if(Math.abs(d) < eps) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!!!!");
				return c;
			} //end if
		} //end for
		
		System.out.println("Max iterations reached without convergence...");
		return c;
	} //end Hybrid method

}
