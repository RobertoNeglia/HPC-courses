#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <numbers>
#include <vector>

class Shape
{
public:
  // constructor that initializes the name
  Shape(const std::string name) : name{name} {}
  // pure virtual getter for the area
  virtual double getArea();
  // getter for the shape name (non virtual)
  const std::string
  getName()
  {
    return name;
  }
  // virtual destructor
  virtual ~Shape() = default;

private:
  // member with the name of the shape (const)
  const std::string name;
};

// Implement the classes "Circle" and "Rectangle"
// The constructor must be empty, in the sense the name of the shape
// should not be a user's choice
class Circle : public Shape
{
public:
  Circle(const double &radius) : Shape("Circle"), radius{radius} {}
  double
  getArea() override
  {
    return radius * radius * 3.14;
  };
  virtual ~Circle() override = default;

private:
  double radius;
};

class Rectangle : public Shape
{
public:
  Rectangle(const double &length, const double &height)
    : Shape("Rectangle"), length(length), height(height){};

  virtual double
  getArea() override
  {
    return length * height;
  };

  virtual ~Rectangle() override = default;

private:
  double length;
  double height;
};

int
main()
{
  // Instantiate vector of shapes
  // Add some shapes
  // Loop over shapes and print
  return 0;
}