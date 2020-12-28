import java.lang.Math; 


public class Find{
	//String[] input;
	Double stoppingCriteria = Math.pow(10,-15);
	
	public static void main(String[] args){
		
	}
	
	public static int factorial(int n){
		// TODO: implement this, not urgent
		return 0;
	}
	
	public static Double singleSolve(double xCur, Func f, ItMe iteration){
		// TODO: implement this
		return null;
	}
	
	public static Double[] completeSolve(Func f, ItMe iteration){
		// TODO: implement this
		return null;
	}
	
	public static Double[] plot(Func f){
		// TODO: implement this
		return null;
	}
}


interface Func{
	public Double evaluate(Double x);
	public Func differentiate();
	
}

interface ItMe{
	public Double next(Double x);
	// Iterative method
	// Turn x_n to x_n+1
}

class Poly implements Func{
	// key * x^ val
	public PolyElement[] elements;
	
	
	public Poly(String input){
		// TODO: implement this
	}
	
	public Poly(PolyElement[] input){
		// TODO: implement this
	}
	
	public Double evaluate(Double x){
		// TODO: implement this
		return null;
	}
	
	public Func differentiate(){
		return null;
	}
	
}

class PolyElement{
	public Double coefficient;
	public int exponent;
	
	PolyElement(Double co, int ex){
		this.coefficient = co;
		this.exponent = ex;
	}
	
	PolyElement(String input){
		// TODO: implement this
		
	}
}















