//
//  importantSampling2D.cpp
//  Mimosa
//
//  Created by Ziyuan Qu on 2023/6/28.
//

#include <iostream>
#include <random>
#include <vector>

// 计算累积分布函数
std::vector<double> computeCDF(const std::vector<std::vector<double>>& data)
{
    std::vector<double> cdf;
    double sum = 0.0;
    for (const auto& row : data)
    {
        for (double value : row)
        {
            sum += value;
            cdf.push_back(sum);
        }
    }
    for (double& value : cdf)
    {
        value /= sum;
    }
    return cdf;
}

// 根据累积分布函数进行重要性采样
std::pair<int, int> importanceSampling(const std::vector<std::vector<double>>& data, const std::vector<double>& cdf)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double r = dis(gen);
    auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
    int index = std::distance(cdf.begin(), it);
    int row = index / data[0].size();
    int col = index % data[0].size();
    return {row, col};
}

int main()
{
    // 定义二维数组
    std::vector<std::vector<double>> data = {{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}, {0.7, 0.8, 0.9}};

    // 计算累积分布函数
    std::vector<double> cdf = computeCDF(data);

    // 进行重要性采样
    auto [row, col] = importanceSampling(data, cdf);

    // 输出结果
    std::cout << "Sampled value: " << data[row][col] << " at (" << row << ", " << col << ")" << std::endl;

    return 0;
}
