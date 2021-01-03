import java.lang.Math; 
import java.util.ArrayList;
import java.util.HashMap;

public class Find{
	//String[] input;
	static Double stoppingCriteria = Math.pow(10,-15);
	
	public static void main(String[] args){
		//Poly hello = new Poly("8:6 7:4 -10:0");
		Poly hello = new Poly("1:5 -56:4 1249:3 -13786:2 75348:1 -163880:0");
		
		println(hello.toString());
		
		Double guess = 9.5;
		
		println("Guess at " + guess.toString());
		Double p = 0.0; //singleSolve(guess, hello, new NewtonRaphson());
		Double q = singleSolve(guess, hello, new HalleyMod());
		Double r = singleSolve(guess, hello, new HouseholderMod());
		Double s = 0.0; //singleSolve(guess, hello, new W4NewtonRaphson(0.5));
		Double t = singleSolve(guess, hello, new DecompositionII());
		
		println(  p.toString()+", "
				+ q.toString()+", " 
				+ r.toString()+", " 
				+ s.toString()+", " 
				+ t.toString());
		
		println("Ended");
	}
	
	public static long factorial(int n){ // https://www.baeldung.com/java-calculate-factorial
	    long fact = 1;
		for (int i = 2; i <= n; i++) {
			fact = fact * i;
		}
		return fact;
	}
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration){
		return singleSolve(xCur, f, iteration, 0);
	}
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration, int depth) throws ArithmeticException{
		// TODO: implement HashMap to avoid calculating the same value twice. 
		// find the next one
		Double xNex = iteration.next(xCur, f);
		depth++;
		
		// determine if a better value is needed
		boolean noImprovement = Math.abs(xNex - xCur) < stoppingCriteria;
		boolean isClose = Math.abs(f.evaluate(xNex)) < stoppingCriteria;
		boolean hasDiverged = depth >= Math.pow(10,3);
		if(hasDiverged){
			throw new ArithmeticException("Diverged");
			//println(Integer.toString(depth));
			//return xNex;
		}
		
		if(noImprovement || isClose){ // other works do not specify if this is AND or OR.
			println(Integer.toString(depth)); // this indicates the quality of the method
			return xNex;
		}else{
			return singleSolve(xNex, f, iteration, depth);
		}
	}
	
	public static Double[] completeSolve(Func f, ItMe iteration){
		// TODO: implement this
		return null;
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
	// Turn x_n to x_n+1
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
	//Double x;
	Double fx  = f.evaluate(x); 
	Double fpx = f.differentiate().evaluate(x);
	Double y   = x - ((fx)/(fpx));
	Double fy  = f.evaluate(y);
	Double fpy = f.differentiate().evaluate(y);
	
	return y - ((2*fx*fy*fpy)/(2*fx*fpy*fpy - fpx*fpx*fy + fpx*fpy*fy));
	}
}

class HouseholderMod implements ItMe{ // Noor et al.: Modified Householder iterative method for nonlinear equations
	HouseholderMod(){}
	
	public Double next(Double x, Func f){
		Double y = new NewtonRaphson().next(x, f);
		Double fy = f.evaluate(y); 
		Double fpy = f.differentiate().evaluate(y);
		Double fppy = f.differentiate().differentiate().evaluate(y);
		
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
	//Double x;
	Double fx  = f.evaluate(x); 
	Double fpx = f.differentiate().evaluate(x);
	Double y   = x - ((fx)/(fpx));
	Double fy  = f.evaluate(y);
	Double fpy = f.differentiate().evaluate(y);
		
		return x- fx/fpx -2*(fy/fpx) + (fy*fpy)/(fpx*fpx);
		
	}
	
}

class Poly implements Func{
	public PolyElement[] elements;
	
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
		// TODO: make LaTeX variant
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
}

class PolyElement implements Func{
	public Double coefficient;
	public int exponent;
	
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
