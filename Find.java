import java.lang.Math; 
import static java.lang.Math.pow;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.io.File;

public class Find{
	static Double stoppingCriteria = Math.pow(10,-15); //Math.sqrt(1.11 * Math.pow(10,-16)); 
	static ArrayList<Double> visitedXValues = new ArrayList<Double>(); 
	
	public static void main(String[] args){
		println("Rootfinding from file:");
		solveFile(args[0]);
		println("Ended");
	}
	
	public static Double safeSolve(Double xCur, Func f, ItMe iteration){
		try{
			return singleSolve(xCur, f, iteration);
		}catch(Exception e){
			System.out.println(e);
			return 0.0;
		}
	}
	
	public static Double expectQuadraticConvergence(Func f, Double xStart){
		Func fp = f.differentiate();
		Func fpp = fp.differentiate();
		
		return Math.abs(f.evaluate(xStart) * fpp.evaluate(xStart)) /(fp.evaluate(xStart)*fp.evaluate(xStart));
	}
	
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration, int depth) throws ArithmeticException{
		// find the next approximation
		Double xNex = iteration.next(xCur, f);
		depth++;
		
		if(visitedXValues.contains(xNex) && Math.abs(f.evaluate(xNex)) > Math.sqrt(stoppingCriteria)){
			//println(xNex.toString());
			throw new ArithmeticException("Non-converging cycle"); // , "Non-converging cycle, x = "+xNex.toString() + ", Steps: " + depth
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
			Double rate = convergenceRate(visitedXValues.toArray(new Double[0]));
			int pop = -1;
			if(f instanceof Counter){
				Counter c = (Counter) f;
				//println("diff.+ eval.: "+Integer.toString(c.popCount()));
				pop = c.popCount();
			}
			/*
			println("Algorithm: " + iteration.toString() +
					"\nSteps: " + Integer.toString(depth) +
					"\nEval.+diff.: " + Integer.toString(pop) +
					"\nConvergence rate: "+ rate //+
					//"\nValue: |f(x_n)| = " + Math.abs(f.evaluate(xNex)) + 
					//"\nRoot: x = " +  Double.toString(xNex)
					);
			println("");
			*/
			
			//String layout = "%s & $%s$ & $%s$ & $%f$ \\\\";
			//layout = String.format(layout, iteration.toString(), Integer.toString(depth), Integer.toString(pop), rate);
			String layout = "& $%s$ & $%s$ & $%f$";
			layout = String.format(layout, Integer.toString(depth), Integer.toString(pop), rate);		

			println(layout);
			
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
							new W4NewtonRaphson(0.5),
							new HalleyMod(), 
							new HouseholderMod(),
							//new DecompositionII(), /* TODO: add back into system*/
							new DecompositionIII(),
							new VariantNewtonsMethod(), 
							new improvedHouseholder(),
							new improvedHouseholderNumerical()
						};

		Double[] roots = new Double[methods.length];
		
		for(int i = 0; i<methods.length ; i++){
			roots[i] = safeSolve(xCur, f, methods[i]);
		}
		return roots;
	}
	
	// Sennning, J. 2007 p. 3
	public static Double convergenceRate(Double mTwo, Double mOne, Double pZero, Double pOne){
		// approximates the convergence rate of a method given 4 entries in the sequence.
		Double a = Math.log(Math.abs((pOne-pZero)/(pZero-mOne)));
		Double b = Math.log(Math.abs((pZero-mOne)/(mOne-mTwo)));
		return a/b;
	}
	public static Double convergenceRate(Double[] sequence){
		int len = sequence.length;
		if (len >=4){
			return convergenceRate(sequence[len-4], sequence[len-3], sequence[len-2], sequence[len-1]);
		}else{
			return -1.0;
		}
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
		Func f = new PolyWithCounter(segmentation[0]); 
		//Func f = new Poly(segmentation[0]);
		String quadratic = Double.toString(expectQuadraticConvergence(f, xStart));
		println("f(x) = "+f.toString()+
				"\nInitial value: x_0 = " + Double.toString(xStart) +
				"\nh(x) = " + quadratic + 
				"\n"
				);

		
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
	// Find x_{n+1} given x_n 
}

class NewtonRaphson implements ItMe{
	NewtonRaphson(){}
	
	public Double next(Double x, Func f){
		return x - ((f.evaluate(x))/(f.differentiate().evaluate(x)));
	}

	public String toString(){
		return "Newton-Raphson";
	}
}

class HalleyMod implements ItMe{ // Noor et al.: A new modified Halley method without second derivatives for nonlinear equation
	private NewtonRaphson nr;
	
	HalleyMod(){
		nr = new NewtonRaphson();
	}

	public Double next(Double x, Func f){
		Double fx  =  f.evaluate(x);
		Func fp = f.differentiate();
		Double fpx =  fp.evaluate(x);
		Double y   =  this.nr.next(x, f);
		Double fy  =  f.evaluate(y);
		Double fpy =  fp.evaluate(y);
		return y - ((2*fx*fy*fpy)/(2*fx*fpy*fpy - fpx*fpx*fy + fpx*fpy*fy));
	}
	public String toString(){
		return "Modified Halley";
	}
}

class HouseholderMod implements ItMe{ // Noor et al.: Modified Householder iterative method for nonlinear equations
	private NewtonRaphson nr;
	
	HouseholderMod(){
		nr = new NewtonRaphson();
	}
	public Double next(Double x, Func f){
		Double y = this.nr.next(x, f);
		Double fy = f.evaluate(y); 
		Func fp = f.differentiate();
		Func fpp = fp.differentiate();
		Double fpy = fp.evaluate(y);
		Double fppy = fpp.evaluate(y);
		
		return y - (fy / fpy) - ((fy*fy*fppy)/(2*fpy*fpy*fpy));
	}
	
	public String toString(){
		return "Modified Householder";
	}
}

class W4NewtonRaphson implements ItMe{ // The W4 method: a new multi-dimensional root-finding scheme for nonlinear systems of equations
	private Double damper; // if this is unchanged, it is indentical to Newton-Raphson
	
	W4NewtonRaphson(Double damp){
		this.damper = damp;
	}
	W4NewtonRaphson(){
		this.damper = 1.0;
	}
	
	public Double next(Double x, Func f){
		return x - damper * ((f.evaluate(x))/(f.differentiate().evaluate(x)));
	}
	
	public String toString(){
		return "W4 method";
	}
}

class DecompositionII implements ItMe{ //Chun, C.:Iterative methods improving Newton's method by the decomposition method 
	private NewtonRaphson nr;
	DecompositionII(){
		nr = new NewtonRaphson();
	}
	
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x);
		Func fp = f.differentiate();
		Double fpx = fp.evaluate(x);
		Double y   = this.nr.next(x, f);
		Double fy  = f.evaluate(y);
		Double fpy = fp.evaluate(y);
		
		return x- fx/fpx -2*(fy/fpx) + (fy*fpy)/(fpx*fpx);
	}
	
	public String toString(){
		return "Decomposition Method II";
	}
}

class DecompositionIII implements ItMe{ //Chun, C.:Iterative methods improving Newton's method by the decomposition method 
	private NewtonRaphson nr;
	DecompositionIII(){
		nr = new NewtonRaphson();
	}
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x);
		Func fp = f.differentiate();
		Func fpp = fp.differentiate();
		Double fpx = fp.evaluate(x);
		Double y   = this.nr.next(x, f);
		Double fy  = f.evaluate(y);
		Double fpy = fp.evaluate(y);		
		Double fppy = fpp.evaluate(y);
		
		return x - (fx)/(fpx) - (3.0*fy)/(fpx) + (3.0*fy*fpy)/(fpx*fpx) - (0.5 * (fy * (fy*fppy+2.0*fpy*fpy)))/(fpy*fpy*fpy);
	}
	public String toString(){
		return "Decomposition Method III";
	}
}

class VariantNewtonsMethod implements ItMe{ // Weerakoon, S.: A variant of Newton's method with accelerated third-order convergence
	private NewtonRaphson nr;
	VariantNewtonsMethod(){
		nr = new NewtonRaphson();
	}
	
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x); 
		Func fp = f.differentiate();
		Double fpx = fp.evaluate(x);
		Double y   = this.nr.next(x, f);
		Double fpy = fp.evaluate(y);
		
		return x- (2*fx)/(fpx + fpy);
	}
	public String toString(){
		return "Variant of Newton's";
	}
}

class improvedHouseholder implements ItMe{// Nazeer, W.: A new Householder method free from second derivatives...
	// third order convergence, but does not need to find f''(x)
	private NewtonRaphson nr;
	improvedHouseholder(){
		nr = new NewtonRaphson();
	}	
	
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x); 
		Func fp = f.differentiate();
		Double fpx = fp.evaluate(x);
		Double y   = this.nr.next(x, f);
		Double fy  = f.evaluate(y);
		Double fpy = fp.evaluate(y);
		
		return y - ((fy)/(fpy)) * (1- (fpy*fpx*fy - fpx*fpx*fx)/(2*fpy*fpy*fx));
	}
	
	public String toString(){
		return "New Householder method";
	}
}

class improvedHouseholderNumerical implements ItMe{
	private Double numericalDerivative(Double x, Func f){
		// Calculates f'(x) for a specific x, avoiding analytical derivatives and its rules.
		// Is significantly less accurate than an analytical solution, but it solves the problem 
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
	
	public String toString(){
		return "New Householder method with numerical derivatives";
	}
}


class Poly implements Func{
	protected PolyElement[] elements;
	
	public Poly(String input){
		input = input.replace("^x", ":");
		input = input.replace("x", ":1");
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

interface Counter{
	public int popCount();
}

class PolyWithCounter extends Poly implements Func, Counter{
	private int count;

	public PolyWithCounter(String input){
		super(input);
		this.count = 0;
	}
	
	public PolyWithCounter(PolyElement[] input){
		super(input);
		this.count = 0;
	}
	
	public int popCount(){
		int tempCount = count;
		this.count = 0;
		return tempCount;
	}
	
	public Double evaluate(Double x){
		this.count++;
		return super.evaluate(x);
	}
	
	public Func differentiate(){
		this.count++;
		return super.differentiate();
	}
}









