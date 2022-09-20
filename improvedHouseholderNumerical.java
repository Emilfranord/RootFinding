import static java.lang.Math.pow;

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