#pragma once

#include <utility>
#include <memory>
#include <vector>

template <typename T>
using vector_ptr = std::shared_ptr<std::vector<T>>;

template <typename T, typename ...Args>
vector_ptr<T> make_vector_ptr(Args &&...args)
{
    return std::make_shared<std::vector<T>>(
        std::forward<Args>(args)...);
}

template <typename Point>
using Polygon = std::tuple<
    vector_ptr<Point>,
    std::vector<unsigned>>;

template <typename Point>
using PolygonPair = std::tuple<
    vector_ptr<Point>,
    std::vector<unsigned>,
    std::vector<unsigned>>;

template <typename Point>
using PolygonMesh = std::tuple<
    vector_ptr<Point>,
    vector_ptr<std::vector<unsigned>>>;

template <typename T, typename Info>
class Decorated: public T
{
    Info info;

public:
    using T::T;

    Decorated(T const &t, Info const &info_):
        T(t), info(info_) {}

    Decorated(Info const &info_):
        info(info_) {}

    Info const &get_info() const
    {
        return info;
    }
};

template <typename Point, typename Info>
using DecoratedMesh = std::tuple<
    vector_ptr<Point>,
    vector_ptr<Decorated<std::vector<unsigned>, Info>>>;
