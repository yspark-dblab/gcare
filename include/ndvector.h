#ifndef NDVECTOR_H_
#define NDVECTOR_H_

#include <algorithm>
#include <vector>

using std::vector;

class Vector1D {
 public:
  size_t size_;
  int* array_;
  Vector1D();
  int operator[](int);
  size_t size(void);
  int* begin(void);
  int* end(void);

  void swap(Vector1D &other) {
    size_t size = size_;
    int* array = array_;
    size_ = other.size_;
    array_ = other.array_;
    other.size_ = size;
    other.array_ = array;
  }
};

class Vector2D {
 public:
  size_t size_;
  int* array_;
  Vector2D();
  class Vector1D {
   public:
    size_t size_;
    int* array_;
    Vector1D(size_t, int*);
    int operator[](int);
    size_t size(void);
    int* begin(void);
    int* end(void);
  };
  Vector1D operator[](int);
  size_t size(void);
  int* begin(void);
  int* end(void);

  void swap(Vector2D &other) {
    size_t size = size_;
    int* array = array_;
    size_ = other.size_;
    array_ = other.array_;
    other.size_ = size;
    other.array_ = array;
  }

};

class Vector3D {
 public:
  size_t size_;
  int* array_;
  Vector3D();
  Vector3D(size_t, int*);
  class Vector2D {
   public:
    size_t size_;
    int* array_;
    Vector2D(size_t, int*);
    class Vector1D {
     public:
      size_t size_;
      int* array_;
      Vector1D(size_t, int*);
      int operator[](int);
      size_t size(void);
      int* begin(void);
      int* end(void);
    };
    Vector1D operator[](int);
    size_t size(void);
    int* begin(void);
    int* end(void);
  };
  Vector2D operator[](int);
  size_t size(void);
  int* begin(void);
  int* end(void);
};

class Vector4D {
 public:
  size_t size_;
  int* array_;
  Vector4D(void);
  Vector4D(size_t, int*);
  class Vector3D {
   public:
    size_t size_;
    int* array_;
    Vector3D(void);
    Vector3D(size_t, int*);
    class Vector2D {
     public:
      size_t size_;
      int* array_;
      Vector2D(size_t, int*);
      class Vector1D {
       public:
        size_t size_;
        int* array_;
        Vector1D(size_t, int*);
        int operator[](int);
        size_t size(void);
        int* begin(void);
        int* end(void);
      };
      Vector1D operator[](int);
      size_t size(void);
      int* begin(void);
      int* end(void);
    };
    Vector2D operator[](int);
    size_t size(void);
    int* begin(void);
    int* end(void);
  };
  Vector3D operator[](int);
  size_t size(void);
  int* begin(void);
  int* end(void);
};

class Vector2DPair {
 public:
  size_t size_;
  int* array_;
  Vector2DPair();
  class Vector1DPair {
   public:
    size_t size_;
    int* array_;
    Vector1DPair(size_t, int*);
    class Pair {
     public:
      int* array_;
      int first, second;
      Pair(int*);
    };
    Pair operator[](int);
    size_t size(void);
    int* begin(void);
    int* end(void);
  };
  Vector1DPair operator[](int);
  size_t size(void);
  int* begin(void);
  int* end(void);
};
#endif

