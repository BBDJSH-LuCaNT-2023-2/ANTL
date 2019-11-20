#include <iostream>
#include <functional>
#include <cstdint>

class foo;

template<>
struct std::hash<foo> {
  public:
    size_t operator()(const foo & f) const;
};

class foo
{
  friend struct std::hash<foo>;

  public:
  foo(unsigned int x, unsigned int y) : x(x), y(y) {}

  private:
  unsigned int x, y;
};

size_t std::hash<foo>::operator()(const foo & f) const
{
  unsigned long z = ((unsigned long)f.x << 32) | (unsigned long)f.y;
  std::hash<unsigned long> h;
  return h(z);
}

int main(int argc, char * argv[])
{

  foo f(1,1);
  std::hash<foo> h;
  std::hash<unsigned long> hl;

  std::cout << h(f) << std::endl;
  for (unsigned long i = 0; i < 10; i++)
    std::cout << hl(i) << std::endl;

  return 0;
}
