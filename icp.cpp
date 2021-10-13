#include <iostream>
using namespace std;
#include "Eigen/Core"
#include "Eigen/SVD"
using namespace Eigen;
#include <cmath>
#include <fstream>
#include <vector>
#define rad_to_deg(rad) (((rad) / 2 / M_PI) * 360)
using namespace std;

int main(int argc, char *argv[])
{

    //配列に入れ直す。
    vector<pair<double, double> > data;
    vector<pair<double, double> > map;

    // MAPデータを読み込み
    FILE *MAP;
    char MAP_file[] = "/home/miura-t/catkin_ws/sample.txt";

    if ((MAP = fopen(MAP_file, "r")) == NULL)
    {
        fprintf(stderr, "%s\n", "error: can't read file.");
        return EXIT_FAILURE;
    }
    double x, y, z;
    while (fscanf(MAP, "%lf %lf", &x, &y) != EOF)
    {
        pair<double, double> p;
        p.first = x;
        p.second = y;
        map.push_back(p);
    }
    fclose(MAP);
    // DATAを読み込み
    FILE *DATA;
    char DATA_file[] = "/home/miura-t/catkin_ws/sample_after.txt";
    if ((DATA = fopen(DATA_file, "r")) == NULL)
    {
        fprintf(stderr, "%s\n", "error: can't read file.");
        return EXIT_FAILURE;
    }
    while (fscanf(DATA, "%lf %lf", &x, &y) != EOF)
    {
        pair<double, double> p;
        p.first = x;
        p.second = y;
        data.push_back(p);
    }
    fclose(DATA);

    //初期化
    MatrixXd data_matrix(2, data.size());
    MatrixXd map_matrix(2, data.size());
    for (int i = 0; i < data.size(); i++)
    {
        data_matrix(0, i) = data[i].first;
        data_matrix(1, i) = data[i].second;
    }

    double pre_err = 10;
    double err;
    int total = 0;
    double translation_x = 0; //合計並進移動量
    double translation_y = 0;
    double spin_rad = 0;                       //合計回転量
    MatrixXd data_move_matrix(2, data.size()); //初期位置で動かした後の点群
    MatrixXd spin_matrix(2, 2);                //特異値分解後に取得するための回転行列
    MatrixXd total_spin_matrix(2, 2);          //初期位置で合わせるための行列
    //角度合わせ
    //初期座標位置合わせ
    total_spin_matrix(0, 0) = cos(spin_rad);
    total_spin_matrix(1, 1) = cos(spin_rad);
    total_spin_matrix(0, 1) = -sin(spin_rad);
    total_spin_matrix(1, 0) = sin(spin_rad);
    for (int i = 0; i < data.size(); i++)
    {
        data_move_matrix(0, i) = translation_x + data_matrix(0, i);
        data_move_matrix(1, i) = translation_y + data_matrix(1, i);
    }
    data_move_matrix = total_spin_matrix * data_move_matrix;

    ///////////////////////////////////////////////////////
    //これ以降繰り返す
    while (total == 0)
    {

        //データ点から最も近い地図データの点を探す。
        double data_x, data_y;
        for (int i = 0; i < data.size(); i++)
        {
            data_x = data_move_matrix(0, i);
            data_y = data_move_matrix(1, i);
            double min_dis, dis;
            min_dis = 100;
            for (int j = 0; j < map.size(); j++)
            {
                dis = (data_x - map[j].first) * (data_x - map[j].first) + (data_y - map[j].second) * (data_y - map[j].second);
                if (dis < min_dis)
                {
                    min_dis = dis;
                    //MAPから取得
                    map_matrix(0, i) = map[j].first;
                    map_matrix(1, i) = map[j].second;
                }
            }
        }
        // 重心を計算する
        //dataの重心
        double ave_data_x = 0;
        double ave_data_y = 0;
        //mapの重心
        double ave_map_x = 0;
        double ave_map_y = 0;
        for (int i = 0; i < data.size(); i++)
        {
            ave_data_x = ave_data_x + data_move_matrix(0, i);
            ave_data_y = ave_data_y + data_move_matrix(1, i);
            ave_map_x = ave_map_x + map_matrix(0, i);
            ave_map_y = ave_map_y + map_matrix(1, i);
        }

        ave_data_x = ave_data_x / data.size();
        ave_data_y = ave_data_y / data.size();
        ave_map_x = ave_map_x / data.size();
        ave_map_y = ave_map_y / data.size();

        //共分散行列の計算
        MatrixXd covariance_data(2, data.size());
        MatrixXd covariance_map(2, data.size());
        for (int i = 0; i < data.size(); i++)
        {
            covariance_data(0, i) = data_move_matrix(0, i) - ave_data_x;
            covariance_data(1, i) = data_move_matrix(1, i) - ave_data_y;
            covariance_map(0, i) = map_matrix(0, i) - ave_map_x;
            covariance_map(1, i) = map_matrix(1, i) - ave_map_y;
        }
        //共分散行列を計算
        MatrixXd covariance(2, 2);
        covariance = covariance_map * covariance_data.transpose() / data.size();

        //共分散行列を特異値分解
        JacobiSVD<MatrixXd> SVD(covariance, ComputeFullU | ComputeFullV);

        //////////////////////////////////////////////////
        MatrixXd s_mat(2, 2);
        MatrixXd M2(2, 2);

        s_mat(0, 0) = 1;
        s_mat(0, 1) = 0;
        s_mat(1, 0) = 0;
        M2 = SVD.matrixU() * SVD.matrixV().transpose();
        s_mat(1, 1) = M2(0, 0) * M2(1, 1) - M2(0, 1) * M2(1, 0);

        spin_matrix = SVD.matrixU() * s_mat * SVD.matrixV().transpose();
        //////////////////////////////////////////////////////
        //並進ベクトルを計算
        MatrixXd ave_map(2, 1);
        MatrixXd ave_data(2, 1);
        ave_map(0, 0) = ave_map_x;
        ave_map(1, 0) = ave_map_y;
        ave_data(0, 0) = ave_data_x;
        ave_data(1, 0) = ave_data_y;

        MatrixXd translation_matrix(2, 1);
        translation_matrix = ave_map - spin_matrix * ave_data;

        cout << "err" << err << endl;

        //誤差の計算
        err = 0;
        double err2 = 0;
        for (int i = 0; i < data.size(); i++)
        {
            err2 = (data_move_matrix(0, i) - map_matrix(0, i)) * (data_move_matrix(0, i) - map_matrix(0, i)) + (data_move_matrix(1, i) - map_matrix(1, i)) * (data_move_matrix(1, i) - map_matrix(1, i));
            err = err + sqrt(err2);
        }
        err = err / data.size();

        //合計移動量、回転角度の更新
        translation_x = translation_x + translation_matrix(0, 0);
        translation_y = translation_y + translation_matrix(1, 0);

        spin_rad = spin_rad + asin(spin_matrix(1, 0));

        std::cout << " 111 " << std::endl;
        if (abs(pre_err - err) < 0.0001)
        {
            std::cout << "END" << std::endl;
            total = 1;
        }
        pre_err = err;
    }

    std::ofstream ofile("/home/miura-t/catkin_ws/icp_result.txt", std::ios::trunc);
    for (int i = 0; i < data.size(); i++)
    {
        ofile << data_move_matrix(0, i) << " " << data_move_matrix(1, i) << " " << 0 << endl;
    }
    ofile.close();

    std::cout << "[x(mm),y(mm),theta(°)]=[" << translation_x * 1000 << "," << translation_y * 1000 << "," << rad_to_deg(spin_rad) << "]"
              << " , "
              << "score[mm]:" << err * 1000
              << std::endl;
}
