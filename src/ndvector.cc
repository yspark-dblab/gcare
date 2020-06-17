#include "../include/ndvector.h"

#include <algorithm>
#include <vector>

using std::vector;
using std::pair;

Vector2D::Vector2D() : size_(0), array_(nullptr) { }
Vector2D::Vector1D::Vector1D(size_t size, int* array) : size_(size), array_(array) { }

Vector1D::Vector1D() : size_(0), array_(nullptr) { }

int Vector1D::operator[](int idx) {
  return array_[idx];
}

size_t Vector1D::size(void) {
  return size_;
}

int* Vector1D::begin(void) {
  return array_;
}

int* Vector1D::end(void) {
  return array_ + size_;
}

int Vector2D::Vector1D::operator[](int idx) {
  return array_[idx];
}

size_t Vector2D::Vector1D::size(void) {
  return size_;
}

int* Vector2D::Vector1D::begin(void) {
  return array_;
}

int* Vector2D::Vector1D::end(void) {
  return array_ + size_;
}

Vector2D::Vector1D Vector2D::operator[](int idx) {
  return Vector2D::Vector1D(array_[idx + 1] - array_[idx] + 1, array_ + idx + array_[idx]);
}

size_t Vector2D::size(void) {
  return size_;
}

int* Vector2D::begin(void) {
  return array_;
}

int* Vector2D::end(void) {
  return array_ + size_;
}

Vector3D::Vector3D() { }
Vector3D::Vector2D::Vector2D(size_t size, int* array) : size_(size), array_(array) { }
Vector3D::Vector2D::Vector1D::Vector1D(size_t size, int* array) : size_(size), array_(array) { }

int Vector3D::Vector2D::Vector1D::operator[](int idx) {
  return array_[idx];
}

size_t Vector3D::Vector2D::Vector1D::size(void) {
  return size_;
}

int* Vector3D::Vector2D::Vector1D::begin(void) {
  return array_;
}

int* Vector3D::Vector2D::Vector1D::end(void) {
  return array_ + size_;
}

Vector3D::Vector2D::Vector1D Vector3D::Vector2D::operator[](int idx) {
  return Vector3D::Vector2D::Vector1D(array_[idx + 1] - array_[idx] + 1, array_ + idx + array_[idx]);
}

size_t Vector3D::Vector2D::size(void) {
  return size_;
}

int* Vector3D::Vector2D::begin(void) {
  return array_;
}

int* Vector3D::Vector2D::end(void) {
  return array_ + size_;
}

Vector3D::Vector2D Vector3D::operator[](int idx) {
  return Vector3D::Vector2D(array_[idx + 1] - array_[idx] + 1, array_ + idx + array_[idx]);
}

size_t Vector3D::size(void) {
  return size_;
}

int* Vector3D::begin(void) {
  return array_;
}

int* Vector3D::end(void) {
  return array_ + size_;
}

Vector4D::Vector4D(void) {
	size_ = 0;
	array_ = nullptr;
}

Vector4D::Vector4D(size_t size, int* array) : size_(size), array_(array) { }
Vector4D::Vector3D::Vector3D(size_t size, int* array) : size_(size), array_(array) { }
Vector4D::Vector3D::Vector2D::Vector2D(size_t size, int* array) : size_(size), array_(array) { }
Vector4D::Vector3D::Vector2D::Vector1D::Vector1D(size_t size, int* array) : size_(size), array_(array) { }

int Vector4D::Vector3D::Vector2D::Vector1D::operator[](int idx) {
	return array_[idx];
}

size_t Vector4D::Vector3D::Vector2D::Vector1D::size(void) {
	return size_;
}

int* Vector4D::Vector3D::Vector2D::Vector1D::begin(void) {
	return array_;
}

int* Vector4D::Vector3D::Vector2D::Vector1D::end(void) {
	return array_ + size_;
}

Vector4D::Vector3D::Vector2D::Vector1D Vector4D::Vector3D::Vector2D::operator[](int idx) {
	return Vector4D::Vector3D::Vector2D::Vector1D(array_[idx + 1] - array_[idx] + 1, array_ + idx + array_[idx]);
}

size_t Vector4D::Vector3D::Vector2D::size(void) {
	return size_;
}

int* Vector4D::Vector3D::Vector2D::begin(void) {
	return array_;
}

int* Vector4D::Vector3D::Vector2D::end(void) {
	return array_ + size_;
}

Vector4D::Vector3D::Vector2D Vector4D::Vector3D::operator[](int idx) {
	return Vector4D::Vector3D::Vector2D(array_[idx + 1] - array_[idx] + 1, array_ + idx + array_[idx]);
}

size_t Vector4D::Vector3D::size(void) {
	return size_;
}

int* Vector4D::Vector3D::begin(void) {
	return array_;
}

int* Vector4D::Vector3D::end(void) {
	return array_ + size_;
}

Vector4D::Vector3D Vector4D::operator[](int idx) {
	return Vector4D::Vector3D(array_[idx + 1] - array_[idx] + 1, array_ + idx + array_[idx]);
}

size_t Vector4D::size(void) {
	return size_;
}

int* Vector4D::begin(void) {
	return array_;
}

int* Vector4D::end(void) {
	return array_ + size_;
}





Vector2DPair::Vector2DPair() : array_(nullptr) { }

Vector2DPair::Vector1DPair::Vector1DPair(size_t size, int* array) : size_(size), array_(array) { }

Vector2DPair::Vector1DPair::Pair::Pair(int* array) : array_(array) { 
  first = array[0];
  second = array[1];
}

Vector2DPair::Vector1DPair::Pair Vector2DPair::Vector1DPair::operator[](int idx) {
  return Pair(array_ + (idx << 1));
}

size_t Vector2DPair::Vector1DPair::size(void) {
  return size_;
}

int* Vector2DPair::Vector1DPair::begin(void) {
  return array_;
}

int* Vector2DPair::Vector1DPair::end(void) {
  return array_ + (size_ << 1);
}

Vector2DPair::Vector1DPair Vector2DPair::operator[](int idx) {
  return Vector2DPair::Vector1DPair((array_[idx + 1] - array_[idx] + 1) >> 1, array_ + idx + array_[idx]);
}

size_t Vector2DPair::size(void) {
  return size_;
}

int* Vector2DPair::begin(void) {
  return array_;
}

int* Vector2DPair::end(void) {
  return array_ + array_[0] - 1;
}

