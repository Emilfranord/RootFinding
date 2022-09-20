import java.util.*;
import java.io.File;

public class Find{
	static Double stoppingCriteria = Math.pow(10,-15); //Math.sqrt(1.11 * Math.pow(10,-16)); 
	static ArrayList<Double> visitedXValues = new ArrayList<Double>(); 
	
	public static void main(String[] args){
		solveManyStandardInput();
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
			visitedXValues.clear();
			throw new ArithmeticException("Non-converging cycle");
		}
		
		visitedXValues.add(xNex);
		
		// determine if a better value is needed
		boolean noImprovement = Math.abs(xNex - xCur) < stoppingCriteria;
		boolean isClose = Math.abs(f.evaluate(xNex)) < stoppingCriteria;
		boolean hasDiverged = depth >= Math.pow(10,4);
		if(hasDiverged){
			visitedXValues.clear();
			throw new ArithmeticException("Diverged");
		}
		
		if(noImprovement || isClose){
			Double rate = convergenceRate(visitedXValues.toArray(new Double[0]));
			int pop = -1;
			
			println("Algorithm: " + iteration.toString() +
					"\nSteps: " + Integer.toString(depth) +
					"\nEval.+diff.: " + Integer.toString(pop) +
					"\nConvergence rate: "+ rate +
					"\nValue: |f(x_n)| = " + Math.abs(f.evaluate(xNex)) + 
					"\nRoot: x = " +  Double.toString(xNex)
					);
			println("");
			
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
							new DecompositionII(), /* removed in favor of DecompositionIII*/
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
		Scanner sr = null;
		try{
			sr = new Scanner(new File(filePath)).useDelimiter("\n");
		} catch(Exception e){
			return null;
		}finally{
			sr.close();
		}
		
		String line = sr.next();
		String[] segmentation = line.split(",");
		
		Double xStart = Double.parseDouble(segmentation[1]);
		Func f = new Poly(segmentation[0]);
		String quadratic = Double.toString(expectQuadraticConvergence(f, xStart));
		println("f(x) = "+f.toString()+
				"\nInitial value: x_0 = " + Double.toString(xStart) +
				"$\nh(x) = " + quadratic + 
				"$\n"
				);
				
		return solveManyMethods(xStart, f);
	}

	public static Collection<Double[]> solveManyStandardInput(){
		Scanner sr = new Scanner(System.in);
		Collection<Double[]> stack = new Stack<>();
		while(sr.hasNextLine()){
			stack.add(solveSingleLine(sr.nextLine())) ;
		}
		sr.close();
		return stack;
	}

	public static Double[] solveSingleLine(String line){
		String[] segmentation = line.split(",");

		Double xStart = Double.parseDouble(segmentation[1]);
		Func f = new Poly(segmentation[0]);
		String quadratic = Double.toString(expectQuadraticConvergence(f, xStart));
		println("f(x) = "+f.toString()+
				"\nInitial value: x_0 = " + Double.toString(xStart) +
				"$\nh(x) = " + quadratic + 
				"$\n"
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