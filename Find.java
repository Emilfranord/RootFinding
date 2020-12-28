import java.lang.Math; 
import java.util.ArrayList;

public class Find{
	//String[] input;
	static Double stoppingCriteria = Math.pow(10,-15);
	
	public static void main(String[] args){
		PolyElement[] hi = {new PolyElement("7:1"), new PolyElement("4:3")};
		Poly hello = new Poly(hi);
		
		
		println("Ended");
	}
	
	public static int factorial(int n){
		// TODO: implement this, not urgent
		return 0;
	}
	
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration){
		return singleSolve(xCur, f, iteration, 0);
	}
	
	public static Double singleSolve(Double xCur, Func f, ItMe iteration, int depth){
		// find the next one
		Double xNex = iteration.next(xCur, f);
		depth++;
		
		// determine if it is good enough
	
		boolean noImprovement = Math.abs(xNex - xCur) < stoppingCriteria;
		boolean isClose = Math.abs(f.evaluate(xNex)) < stoppingCriteria;
		boolean diverged = depth >= Math.pow(10,4);
		if(diverged){
			println(depth.toString());
			return null;
		}
		
		if(noImprovement && isClose){
			println(depth.toString()); // this indicates the quality of the method
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
}

interface Func{
	public Double evaluate(Double x);
	public Func differentiate();
}

interface ItMe{
	public Double next(Double x, Func f);
	// Iterative method
	// Turn x_n to x_n+1
}


class 


class Poly implements Func{
	public PolyElement[] elements;
	
	public Poly(String input){
		// TODO: implement this
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
