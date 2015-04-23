#ifndef RAND_VAR_HPP
#define RAND_VAR_HPP


/*************************************************************************************
	* STRUCT Rand_Var
	*
	* Constraints on type 'Dist':
	*						- typedef 'result_type'
	*						- static constexpr 'dim_alea' : the number of quasi/pseudo
	*						  random numbers needed to yield
	*						- operator(const Array<Dist::dim_alea>&) yielding 'result_type'
	*
	* Constraint on type 'Generator':
	*						- operator() yielding Array<dim>
	*
	*************************************************************************************/
template<typename Dist, typename Generator>
struct Rand_Var{

				typedef typename Dist::result_type result_type;

				Rand_Var(const Dist& dist, const Generator& gen):
								dist_(dist), gen_(gen), value_(dist_(gen_())) {}

				result_type operator()(){ return value_ = dist_(gen_()); }
				result_type current() const { return value_; }

protected:
				Dist dist_;
				Generator gen_;
				result_type value_;
};



/*************************************************************************************
	* FUNC make_rvar<Dist,Generator>
	*
	* Bind a generator (of real in a unit hypercube) to a distribution.
	*
	*************************************************************************************/
template<typename Dist, typename Generator>
Rand_Var<Dist,Generator>
make_rvar(const Dist& dist, const Generator& gen)
{
				return Rand_Var<Dist,Generator>(dist,gen);
}




























#endif // RAND_VAR_HPP
