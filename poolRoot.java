import java.util.*;
import java.io.*;

public class polRoot {

	static int defaultIter = 10000;
	static float[] solutionArray = new float[3];
	static float a, b;
	static String mode;

	public static void main(String[] args) {
		String[] delim;
		String fileName = "";
		String solutionName = "";
		//> polRoot [-newt, -sec] [-maxIt n] initP [initP2] polyFileName
		float initP1 = 0;
        float initP2 = 0;
        boolean hasInitP = false;
        mode = "bisection";
		for(int i = 0; i < args.length; i++) {
			if(args[i].equalsIgnoreCase("-newt")) {
				System.out.println("Newton Method");
				mode = "newton";
				continue;
			}
			else if(args[i].equalsIgnoreCase("-sec")) {
				System.out.println("Secant Method");
				mode = "secant";
				continue;
			}
			else if(args[i].equalsIgnoreCase("-hybrid")) {
				System.out.println("Hybrid Method");
				mode = "hybrid";
				continue;
			}
			if (args[i].equalsIgnoreCase("-maxIt")) {
                defaultIter = Integer.parseInt(args[i+1]);
                continue;
            }
			if (!args[i].endsWith(".pol")) {
                if (!hasInitP) {
                    initP1 = Float.parseFloat(args[i]);
                    hasInitP = true;
                } else {
                    initP2 = Float.parseFloat(args[i]);
                }

                continue;
            }
			if(args[i].endsWith(".pol")) {
				a = initP1;
				b = initP2;
				fileName = args[i];
				delim = fileName.split("\\.");
				solutionName = delim[0] + "." + "sol";
			}
		}
		if(mode.equalsIgnoreCase("bisection"))
			System.out.println("Bisection Method");
		readAndWrite(fileName,solutionName);
		System.exit(0);
	}
	public static void readAndWrite(String fileName,String solutionName) {
		try {
			BufferedReader bufferedReader = new BufferedReader(new FileReader(System.getProperty("user.dir")+"/"+fileName));

			//Read first line to determine polynomial degree
			int degree = Integer.parseInt(bufferedReader.readLine());
			float[] poly = new float[degree+1];
			//Read in coefficients
			String line = bufferedReader.readLine();
			String[] coeff = line.split("\\s+");
			for(int i = 0; i <= degree; i++) {
				poly[i] = Float.parseFloat(coeff[i]);
			}
		/*
		 * Execute specific methods
		 */
			switch(mode) {
			case "bisection":
				Bisection(poly, a, b, defaultIter, 0);
				break;
			case "newton":
				Newton(poly, DeriveFunction(poly), a, defaultIter, 0, 0);
				break;
			case "secant":
				Secant(poly, a, b, defaultIter, 0);
				break;
			case "hybrid":
				Hybrid(poly, a, b, defaultIter, 0, 0);
				break;
			}
			String answer = Arrays.toString(solutionArray);
			answer = answer.replaceAll("\\[", "").replaceAll("\\]","").replace("\\,", " ");
			writeEntries(answer,solutionName);
			bufferedReader.close();
		}//end try
		catch(FileNotFoundException e) {
		}
		catch(IOException e) {
		}//end catch
	}//end readEntries
	public static void writeEntries(String answer, String solutionFile) {
		//Writer
		try {
			BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(solutionFile));
			bufferedWriter.write(answer);
            bufferedWriter.close();
		}
		catch(IOException e) {
		}
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
			if( Math.abs(error) == eps || fc == 0) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!");
				solutionArray[0] = c;
				solutionArray[1] = iterCounter;
				solutionArray[2] = 1;
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
		solutionArray[0] = c;
		solutionArray[1] = iterCounter;
		solutionArray[2] = -1;
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
			if(Math.abs(d) == eps) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!");
				solutionArray[0] = x;
				solutionArray[1] = iterCounter;
				solutionArray[2] = 1;
				return x;
			} //end if
		} //end for
		System.out.println("Max iterations reached without convergence...");
		solutionArray[0] = x;
		solutionArray[1] = iterCounter;
		solutionArray[2] = -1;
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
			if(Math.abs(d) == eps) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!");
				solutionArray[0] = a;
				solutionArray[1] = iterCounter;
				solutionArray[2] = 1;
				return a;
			} //end if
			a = a - d;
			fa = Function(f, a);
		} //end for
		System.out.println("Maximum number of itertions reached!!");
		solutionArray[0] = a;
		solutionArray[1] = iterCounter;
		solutionArray[2] = -1;
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
			if( Math.abs(error) == eps || fc == 0) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!!");
				solutionArray[0] = c;
				solutionArray[1] = iterCounter;
				solutionArray[2] = 1;
				return c;
			} //end if
			else if(Math.abs(error) >= 0.01) {
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
			//Stop once reached max iteration
			if(iterCounter >= maxIter)
				break;
			iterCounter++;
			float fd = Function(DeriveFunction(f), c);
			if ( Math.abs(fd) < delta) {
				//Continue with bisection method instead by breaking out Newton's for loop
				break;
			} //end if
			float d = (float)(fc / fd);
			//New c
			c = c - d;
			fc = Function(f, c);
			if(Math.abs(d) == eps) {
				System.out.println("Algorithm has converged after " + iterCounter + " iterations!!!");
				solutionArray[0] = c;
				solutionArray[1] = iterCounter;
				solutionArray[2] = 1;
				return c;
			} //end if
		} //end for
		System.out.println("Max iterations reached without convergence...");
		solutionArray[0] = c;
		solutionArray[1] = iterCounter;
		solutionArray[2] = -1;
		return c;
	} //end Hybrid method
}
