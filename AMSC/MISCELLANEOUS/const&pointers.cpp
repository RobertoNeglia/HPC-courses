#include <iostream>

int main() {
  int a = 10;
  int b = 20;

  int const *w = &a;
  // cant do this (change the pointed value)
  // *w = a;
  // but i can do this (change the pointer value)
  w = &b;
  // -> the pointed value cant be changed with this pointer, but the pointer can
  // be changed

  const int *const x = &a;
  // cant do this (change the pointed value)
  // *x = 11;
  // and cant do this either (change pointer value)
  // x = &b;
  // -> the pointer and the pointed value cant both be changed

  const int *y = &a;
  // cant do this
  // *y = 11;
  // but i can do this (change the pointer value)
  y = &b;
  // -> the pointed value cant be changed with this pointer, but the pointer can
  // be changed
  // ----> w and y are the same

  int *const z = &a;
  // can do this (change the pointed value)
  *z = 11;
  // but cant do this (change the pointer value)
  // z = &b;
  // -> the pointed value can be changed, but not the pointer to that variable
}