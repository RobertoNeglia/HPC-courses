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
  Shape(const std::string &name) : name{name} {}
  // pure virtual getter for the area
  // = 0 is necessary since it requires the function to be overwritten in a
  // derived class
  virtual double getArea() = 0;
  // getter for the shape name (non virtual)
  // The const after the method name means that the "this" argument (the object
  // itself) is const qualified: this means that this method can't change the
  // state of the object
  const std::string &getName() const { return name; }
  // virtual destructor
  // always make base classes' destructors virtual when they're meant to be
  // manipulated polymorphically.
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
  Circle(const double radius) : Shape("Circle"), radius{radius} {}
  double getArea() override
  {
    return radius * radius * std::numbers::pi_v<double>;
  };
  // = default -> use compiler generated version of function (in this case
  // destructor)
  virtual ~Circle() override = default;

private:
  const double radius;
};

class Rectangle : public Shape
{
public:
  Rectangle(const double length, const double height)
      : Shape("Rectangle"), length(length), height(height){};

  virtual double getArea() override { return length * height; };

  virtual ~Rectangle() override = default;

private:
  const double length;
  const double height;
};

int main()
{
  // Instantiate vector of shapes
  // Add some shapes
  // Loop over shapes and print
  std::vector<std::shared_ptr<Shape>> shapes;
  shapes.push_back(std::make_shared<Circle>(2.0));
  shapes.push_back(std::make_shared<Rectangle>(1.5, 5.43));
  shapes.push_back(std::make_shared<Circle>(123.34));

  for (const auto s : shapes)
  {
    std::cout << "I am a " << s->getName()
              << ", and my area is: " << s->getArea() << std::endl;
  }
  return 0;
}