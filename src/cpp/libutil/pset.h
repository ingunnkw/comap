#ifndef PSETINC
#define PSETINC

#include <vector>
#include <set>

template <class T, class Container = std::set<T> >
class Pset : public Container
{
public:
	typedef typename Container::iterator iterator;
	typedef typename Container::iterator const_iterator;
	typedef typename Container::value_type value_type;
	Pset() {}
	Pset(const std::vector<T> & v)
	{
		for(int i = 0; i < v.size(); i++)
			Container::insert(v[i]);
	}
	Pset(const Pset & p) { for(const_iterator i = p.begin(); i != p.end(); i++) Container::insert(*i); }
	operator std::vector<int>() const
	{
		std::vector<int> res;
		for(const_iterator i = Container::begin(); i != Container::end(); i++)
			res.push_back(*i);
		return res;
	}
	bool contains(int i) const { return Container::find(i) != Container::end(); }
	bool contains(const Pset & p) const
	{
		for(const_iterator i = p.begin(); i != p.end(); i++)
			if(!contains(*i)) return false;
		return true;
	}
	bool subset_of(const Pset & p) const { return p.contains(*this); }
	// Intersection is much faster if you loop over the smallest set.
	Pset intersection(const Pset & p) const
	{
		Pset res;
		if(Container::size() < p.size())
		{
			for(const_iterator i = Container::begin(); i != Container::end(); i++)
				if(p.contains(*i)) res.insert(*i);
		}
		else
		{
			for(const_iterator i = p.begin(); i != p.end(); i++)
				if(contains(*i)) res.insert(*i);
		}
		return res;
	}
	Pset union_with(const Pset & p) const
	{
		Pset res = *this;
		for(const_iterator i = Container::begin(); i != Container::end(); i++)
			res.insert(*i);
		return res;
	}
	Pset disjunction(const Pset & p) const
	{
		Pset res;
		for(const_iterator i = Container::begin(); i != Container::end(); i++)
			if(!p.contains(*i)) res.insert(*i);
		for(const_iterator i = p.begin(); i != p.end(); i++)
			if(!contains(*i)) res.insert(*i);
		return res;
	}
	Pset operator+(const Pset & p) const { return union_with(p); }
	Pset operator|(const Pset & p) const { return union_with(p); }
	Pset operator-(const Pset & p) const
	{
		Pset res = p;
		for(const_iterator i = Container::begin(); i != Container::end(); i++)
			res.erase(*i);
		return res;
	}
	Pset operator^(const Pset & p) const { return disjunction(p); }
	Pset operator&(const Pset & p) const { return intersection(p); }
	Pset & operator+=(const Pset & p) {
		for(const_iterator i = p.begin(); i != p.end(); i++) Container::insert(*i);
		return *this;
	}
	Pset & operator|=(const Pset & p) { return *this += p; }
	Pset & operator-=(const Pset & p) {
		for(const_iterator i = p.begin(); i != p.end(); i++) Container::erase(*i);
		return *this;
	}
	Pset & operator^=(const Pset & p) { return *this = *this ^ p; }
	Pset & operator&=(const Pset & p) { return *this = *this & p; }
	bool operator<(const Pset & p) const { return subset_of(p); }
	bool operator>(const Pset & p) const { return contains(p); }

	Pset operator+(T i) const { Pset res = *this; return res += i; }
	Pset operator|(T i) const { Pset res = *this; return res |= i; }
	Pset operator-(T i) const { Pset res = *this; return res -= i; }
	Pset operator^(T i) const { Pset res = *this; return res ^= i; }
	Pset & operator+=(T i) { Container::insert(i); return *this; }
	Pset & operator|=(T i) { Container::insert(i); return *this; }
	Pset & operator-=(T i) { erase(i);  return *this; }
	Pset & operator^=(T i) { return contains(i) ? *this -= i : *this += i; }
	bool operator<(T i) const { return Container::empty(); }
	bool operator>(T i) const { return contains(i); }
};

template<class T> bool operator<(T i, const Pset<T> & p) { return p > i; }
template<class T> bool operator>(T i, const Pset<T> & p) { return p < i; }

#endif
