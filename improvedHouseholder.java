
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
		
		//return y - ((fy)/(fpy)) * (1- (fpy*fpx*fy - fpx*fpx*fx)/(2*fpy*fpy*fx));
		return y - ((fy)/(fpy)) * (1- (fpy*fpx*fy - fp.evaluate(fp.evaluate(x))*fx)/(2* fp.evaluate(fp.evaluate(y)) *fx));
		// The second method likely more correct, but i has not been used in the tests.
		// It depends on the syntax f'^2(x), and if that is f'(f'(x)) or f'(x)* f'(x).
	}
	
	public String toString(){
		return "New Householder method";
	}
}