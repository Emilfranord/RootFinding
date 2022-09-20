
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