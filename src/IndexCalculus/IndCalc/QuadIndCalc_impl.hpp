#ifndef QUADINDCALC_IMPL_H
#define QUADINDCALC_IMPL_H

template <class T, class R>
void QuadIndCalc<T,R>::compute_cg_fac_base(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  IndCalc<R>::cg_factor_base = std::vector<IOrder>();
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_cg_relations(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  IndCalc<R>::cg_relations = std::vector<IOrder>();
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_cg_mat(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  IndCalc<R>::cg_mat.SetDims(IndCalc<R>::cg_relations.size(), IndCalc<R>::cg_factor_base.size());
};

#endif //QUADCLASSGROUPINDCALC_H
