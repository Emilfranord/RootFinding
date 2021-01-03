import java.lang.Math; 
import java.util.ArrayList;
import java.util.HashMap;

public class Find{
	//String[] input;
	static Double stoppingCriteria = Math.pow(10,-15);
	
	public static void main(String[] args){
		
		Poly hello = new Poly("8:6 7:4 99:0");
		
		Double p = singleSolve(1.5, hello , new NewtonRaphson());
		Double q = singleSolve(1.5, hello , new HalleyMod());
		Double r = singleSolve(1.5, hello , new HouseholderMod());
		
		println(q.toString() +", "+ p.toString() +", " + r.toString());
		println(hello.toString());
		
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
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration, int depth){
		// TODO: implement HashMap to avoid calculating the same value twice. 
		// find the next one
		Double xNex = iteration.next(xCur, f);
		depth++;
		
		// determine if a better value is needed
		boolean noImprovement = Math.abs(xNex - xCur) < stoppingCriteria;
		boolean isClose = Math.abs(f.evaluate(xNex)) < stoppingCriteria;
		boolean hasDiverged = depth >= Math.pow(10,4);
		if(hasDiverged){
			println("Diverged at: "+Integer.toString(depth));
			return xNex;
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

class W4NewtonRaphson implements ItMe{
	private Double damper = 1; // if this is unchanged, it is indentical to NewtonRaphson
	
	W4NewtonRaphson(Double damp){
		if (damp > 1){
			this.damp = 1;
			println("Error, damp should be in the interval [0;1]")
		}
		this.damper = damp;
	}
	W4NewtonRaphson(){
		this.damper = 1;
	}
	
	public Double next(Double x, Func f){
		
		return x - damper * ((f.evaluate(x))/(f.differentiate().evaluate(x)));
		
	}
	
}



class Poly implements Func{
	public PolyElement[] elements;
	
	public Poly(String input){
		// TODO: implement this
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
		// TODO: improve this
		StringBuilder temp = new StringBuilder();
		for(PolyElement q :elements){
			temp.append(q.toString());
			temp.append("+");
		}
		sb.setLength(sb.length() - 1);
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
	
	public Double evaluate(Double x){
		return this.coefficient * Math.pow(x, this.exponent);
	}
	
	public Func differentiate(){
		return new PolyElement(this.coefficient * this.exponent ,this.exponent-1);
	}
	
}
