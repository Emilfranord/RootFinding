
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