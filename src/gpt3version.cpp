#include <iostream>
#include <vector>

using namespace std;

typedef vector<int> Poly;

const int N = 1005;

int n, m;
Poly f, g, r, s, t;

int p;  // Modulus of the finite field

int inv(int a) {
  int b = p - 2, c = 1;
  while (b) {
    if (b & 1) c = 1ll * c * a % p;
    a = 1ll * a * a % p;
    b >>= 1;
  }
  return c;
}

Poly operator- (const Poly &a, const Poly &b) {
  int n = a.size(), m = b.size();
  Poly c(max(n, m));
  for (int i = 0; i < n; i++) {
    c[i] = a[i];
  }
  for (int i = 0; i < m; i++) {
    c[i] = (c[i] - b[i] + p) % p;
  }
  return c;
}

Poly operator* (const Poly &a, const Poly &b) {
  int n = a.size(), m = b.size();
  Poly c(n + m - 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      c[i + j] = (c[i + j] + a[i] * b[j]) % p;
    }
  }
  return c;
}

Poly operator/ (Poly a, Poly b) {
  int n = a.size(), m = b.size();
  Poly c(n - m + 1);
  int d = inv(b[m - 1]);
  for (int i = n - 1; i >= m - 1; i--) {
    c[i - m + 1] = 1ll * d * a[i] % p;
    for (int j = 0; j < m && i - j - 1 >= 0; j++) {
      a[i - j - 1] = (a[i - j - 1] - 1ll * c[i - m + 1] * b[m - 1 - j] % p + p) % p;
    }
  }
  return c;
}

Poly operator% (Poly a, Poly b) {
  int n = a.size(), m = b.size();
  Poly c = a / b * b;
  for (int i = 0; i < m; i++) {
    a[i] = (a[i] - c[i] + p) % p;
  }
  while (a.size() > m - 1) a.pop_back();
  return a;
}

void exgcd(Poly f, Poly g) {
  if (g.size() == 1) {
    r = {inv(g[0])};
    s = {1};
    t = {0};
    return;
  }
  Poly h = f % g;
  exgcd(g, h);
  Poly p = r;
  r = s;
  s = p - s * (f / g);
}

pair<Poly, Poly> extended_euclidean_algorithm(Poly f, Poly g, int field_size) {
  p = field_size;
  exgcd(f, g);
  return {r, s};
}


int main() {
  // Read in the polynomials f and g and the cardinality of the finite field
  cin >> n >> m >> p;
  f.resize(n);
  g.resize(m);
  for (int i = 0; i < n; i++) cin >> f[i];
  for (int i = 0; i < m; i++) cin >> g[i];

  // Compute the pair (U, V) such that f*U + g*V = gcd(f, g)
  pair<Poly, Poly> result = extended_euclidean_algorithm(f, g, p);
  Poly U = result.first;
  Poly V = result.second;

  // Print the results
  cout << "U = ";
  for (int i = 0; i < U.size(); i++) {
    cout << U[i] << " ";
  }
  cout << endl;
  cout << "V = ";
  for (int i = 0; i < V.size(); i++) {
    cout << V[i] << " ";
  }
  cout << endl;

  return 0;
}

