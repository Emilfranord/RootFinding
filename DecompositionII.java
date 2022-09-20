
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