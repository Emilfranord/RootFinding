import java.lang.Math; 
import static java.lang.Math.pow;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.io.File;

// use flag: -Xss8m, if java throws a stackoverflow

public class Find{
	//String[] input;
	static Double stoppingCriteria = Math.pow(10,-15);
	static ArrayList<Double> visitedXValues = new ArrayList<Double>(); 
	
	public static void main(String[] args){
		//Poly hello = new Poly("8:6 7:4 -10:0");
		//Poly hello = new Poly("1:5 -56:4 1249:3 -13786:2 75348:1 -163880:0");
		
		//println(hello.toString());
		//println(hello.toString("latex"));
		
		//Double guess =19.1;
		
		//println("Guess at " + guess.toString());
		//solveManyMethods(guess, hello);
		
		//println("solve from file:");
		//solveFile("test.txt");
		
		println("solve from file:");
		solveFile(args[0]);

		
		println("Ended");
	}
	
	public static long factorial(int n){ // https://www.baeldung.com/java-calculate-factorial
	    long fact = 1;
		for (int i = 2; i <= n; i++) {
			fact = fact * i;
		}
		return fact;
	}
	
	public static Double safeSolve(Double xCur, Func f, ItMe iteration){
		try{
			return singleSolve(xCur, f, iteration);
		}catch(Exception e){
			System.out.println(e);
			return 0.0;
		}
	}
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration, int depth) throws ArithmeticException{
		// find the next one
		Double xNex = iteration.next(xCur, f);
		depth++;
		
		if(visitedXValues.contains(xNex)){
			throw new ArithmeticException("Non-converging cycle");
		}
		
		visitedXValues.add(xNex);
		
		// determine if a better value is needed
		boolean noImprovement = Math.abs(xNex - xCur) < stoppingCriteria;
		boolean isClose = Math.abs(f.evaluate(xNex)) < stoppingCriteria;
		boolean hasDiverged = depth >= Math.pow(10,4);
		if(hasDiverged){
			throw new ArithmeticException("Diverged");
		}
		
		if(noImprovement || isClose){
			println("Root: x = " +  Double.toString(xNex) + ", Steps: " + Integer.toString(depth)); // depth indicates the quality of the method
			visitedXValues.clear();
			return xNex;
		}else{
			return singleSolve(xNex, f, iteration, depth);
		}
	}
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration){
		return singleSolve(xCur, f, iteration, 0);
	}
	
	public static Double[] solveManyMethods(Double xCur, Func f){
		ItMe[] methods = {	new NewtonRaphson(), 
							new HalleyMod(), 
							new HouseholderMod(),
							new W4NewtonRaphson(0.5),
							new DecompositionII(), 
							new VariantNewtonsMethod(), 
							new improvedHouseholder(),
							new improvedHouseholderNumerical()
						};

		Double[] roots = new Double[methods.length];
		
		
		for(int i = 0; i<methods.length ; i++){
			roots[i]= safeSolve(xCur, f, methods[i]);
			
		}
		return roots;
		
	}
	
	// Sennning, J. 2007 p. 3
	public static Double convergenceRate(Double mTwo,Double mOne, Double pZero, Double pOne){
		// approximates the convergence rate of a method given 4 entries in the sequence.
		Double a = Math.log(Math.abs(pOne-pZero)/(pZero-mOne));
		Double b = Math.log(Math.abs(pZero-mOne)/(mOne-mTwo));
		return a/b;
	}
	public static Double convergenceRate(Double[] sequence){
		int len = sequence.length;
		return convergenceRate(sequence[len-4], sequence[len-3], sequence[len-2], sequence[len-1]);
	}
	
	public static Double[] solveFile(String filePath){
		Scanner sr;
		try{
			sr = new Scanner(new File(filePath)).useDelimiter("\n");
		} catch(Exception e){
			return null;
		}
		
		String line = sr.next();
		String[] segmentation = line.split(",");
		
		Double xStart = Double.parseDouble(segmentation[1]);
		Func f = new Poly(segmentation[0]);
		println("The function in the file is: "+f.toString());
		
		return solveManyMethods(xStart, f);
	}


	public static ArrayList<Double> plot(Func f, Double jumpLength, Double min, Double max ){
		ArrayList<Double> temp = new ArrayList<Double>();
		
		for(Double d = min ; d < max;  d+=jumpLength){
			temp.add(f.evaluate(d));
		}
		return temp;
	}
	
	public static void println(String str){
		System.out.println(str);
	}
	
	public static boolean isDoubleRoot(Func f, Double x){
		return Math.abs(f.differentiate().evaluate(x)) < stoppingCriteria;
	}
}

interface Func{
	public Double evaluate(Double x);
	public Func differentiate();
}

interface ItMe{ // Iterative method
	public Double next(Double x, Func f);
	// Find x_n+1 given x_n 
}

class NewtonRaphson implements ItMe{
	NewtonRaphson(){}
	
	public Double next(Double x, Func f){
		return x - ((f.evaluate(x))/(f.differentiate().evaluate(x)));
	}
}

class HalleyMod implements ItMe{ // Noor et al.: A new modified Halley method without second derivatives for nonlinear equation
	HalleyMod(){}

	public Double next(Double x, Func f){
	Double fx  =  f.evaluate(x);
	Func fp = f.differentiate();
	Double fpx =  fp.evaluate(x);
	Double y   =  x - ((fx)/(fpx));
	Double fy  =  f.evaluate(y);
	Double fpy =  fp.evaluate(y);
	
	return y - ((2*fx*fy*fpy)/(2*fx*fpy*fpy - fpx*fpx*fy + fpx*fpy*fy));
	}
}

class HouseholderMod implements ItMe{ // Noor et al.: Modified Householder iterative method for nonlinear equations
	HouseholderMod(){}
	
	public Double next(Double x, Func f){
		Double y = new NewtonRaphson().next(x, f);
		Double fy = f.evaluate(y); 
		Func fp = f.differentiate();
		Func fpp = fp.differentiate();
		Double fpy = fp.evaluate(y);
		Double fppy = fpp.evaluate(y);
		
		return y - (fy / fpy) - ((fy*fy*fppy)/(2*fpy*fpy*fpy));
	}
}

class W4NewtonRaphson implements ItMe{ // The W4 method: a new multi-dimensional root-finding scheme for nonlinear systems of equations
	private Double damper; // if this is unchanged, it is indentical to NewtonRaphson
	
	W4NewtonRaphson(Double damp){
		this.damper = damp;
	}
	W4NewtonRaphson(){
		this.damper = 1.0;
	}
	
	public Double next(Double x, Func f){
		return x - damper * ((f.evaluate(x))/(f.differentiate().evaluate(x)));
	}
}

class DecompositionII implements ItMe{ //Chun, C.:Iterative methods improving Newton's method by the decomposition method 
	public Double next(Double x, Func f){
	Double fx  = f.evaluate(x); 
	Double fpx = f.differentiate().evaluate(x);
	Double y   = x - ((fx)/(fpx));
	Double fy  = f.evaluate(y);
	Double fpy = f.differentiate().evaluate(y);
		
		return x- fx/fpx -2*(fy/fpx) + (fy*fpy)/(fpx*fpx);
	}
}

class VariantNewtonsMethod implements ItMe{ // Weerakoon, S.: A variant of Newton's method with accelerated third-order convergence
	VariantNewtonsMethod(){}
	
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x); 
		Double fpx = f.differentiate().evaluate(x);
		Double y   = x - ((fx)/(fpx)); // called x* in this article.
		
		Double fpy = f.differentiate().evaluate(y);
		
		return x- (2*fx)/(fpx + fpy);
	}
}

class improvedHouseholder implements ItMe{// Nazeer, W.: A new Householder method free from second derivatives...
	// third order convergence, but does not need to find f''(x)
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x); 
		Double fpx = f.differentiate().evaluate(x);
		Double y   = x - ((fx)/(fpx));
		Double fy  = f.evaluate(y);
		Double fpy = f.differentiate().evaluate(y);
		
		return y - ((fy)/(fpy)) * (1- (fpy*fpx*fy - fpx*fpx*fx)/(2*fpy*fpy*fx));
		
	}
}

class improvedHouseholderNumerical implements ItMe{
	private Double numericalDerivative(Double x, Func f){
		// Calculates f'(x) for a specific x, avoiding analytical derivatives and its rules.
		// TODO: Implement this
		double deltaX;

		for (int exponent = 0; exponent< 32; exponent++) {
			deltaX = pow(10.0, -exponent);
			double numericalValue = (f.evaluate(x+deltaX)- f.evaluate(x))/(deltaX);

			if (numericalValue == 0.0) { 
				// Testing has shown that the most accurate approximation is found
				// by letting deltaX be 100 times bigger, than when the fraction is equal to zero
				exponent--;
				exponent--;
				deltaX = pow(10.0, -exponent);
				return (f.evaluate(x+deltaX)- f.evaluate(x))/(deltaX);
			}
		}
		return  0.0;
	}
	
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x); 
		Double fpx = numericalDerivative(x,f);
		Double y   = x - ((fx)/(fpx));
		Double fy  = f.evaluate(y);
		Double fpy = numericalDerivative(y,f);
		
		return y - ((fy)/(fpy)) * (1- (fpy*fpx*fy - fpx*fpx*fx)/(2*fpy*fpy*fx));
	}
}


class Poly implements Func{
	protected PolyElement[] elements;
	
	public Poly(String input){
		String[] splitInp = input.split(" ");
		PolyElement[] temp = new PolyElement[splitInp.length];
		
		int index = 0;
		for (String q : splitInp){
			temp[index] = new PolyElement(q);
			index++;
		}
		this.elements = temp;
		
	}
	
	public Poly(PolyElement[] input){
		this.elements = input;
	}
	
	public Double evaluate(Double x){
		Double sum = 0.0;
		for(PolyElement q : this.elements){
			sum += q.evaluate(x);
		}
		return sum;
	}
	
	public Func differentiate(){
		PolyElement[] future = new PolyElement[elements.length];
		for(int i= 0;i <elements.length;i++){
			future[i] = (PolyElement) elements[i].differentiate();
		} 
		return new Poly(future);
	}
	
	public String toString(){
		StringBuilder temp = new StringBuilder();
		boolean first = true;
		for(PolyElement q :elements){
			if (q.coefficient >=0){
				if (!first){
					temp.append("+");
				}
			}
			temp.append(q.toString());
			first = false;
		}
		return temp.toString();
	}
	
	public String toString(String type){
		if(!type.equals("latex")){
			return this.toString();
		}
		StringBuilder temp = new StringBuilder();
		boolean first = true;
		for(PolyElement q :elements){
			if (q.coefficient >=0){
				if (!first){
					temp.append("+");
				}
			}
			temp.append(q.toString("latex"));
			first = false;
		}
		return temp.toString();
	}
	
}

class PolyElement implements Func{
	protected Double coefficient;
	protected int exponent;
	
	PolyElement(Double co, int ex){
		this.coefficient = co;
		this.exponent = ex;
	}
	
	PolyElement(String input){
		//expect format like "3:6", or "C:E"
		String[] temp = input.split(":", 2);
		this.coefficient = Double.parseDouble(temp[0]);
		this.exponent = Integer.parseInt(temp[1]);
	}
	
	public String toString(){
		if(exponent == 0){return coefficient.toString();}
		if(exponent == 1){return coefficient +"*x";}
		if(coefficient == 0){return "";}
		
		return coefficient +"*x^" + exponent;
	}
	
	public String toString(String type){
		if(!type.equals("latex")){
			return this.toString();
		}
		
		if(exponent == 0){return coefficient.toString();}
		if(exponent == 1){return coefficient +"\\cdot x";}
		if(coefficient == 0){return "";}
		
		return coefficient +"\\cdot x^{" + exponent+"}";
	}
	
	public Double evaluate(Double x){
		return this.coefficient * Math.pow(x, this.exponent);
	}
	
	public Func differentiate(){
		return new PolyElement(this.coefficient * this.exponent ,this.exponent-1);
	}
	
}
