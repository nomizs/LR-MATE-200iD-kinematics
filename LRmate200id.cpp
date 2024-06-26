#include <iostream>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Dense>      //biblioteka umo¿liwiaj¹ca tworzenie oraz wykonywanie dzia³an na macierzach i wektorach
#include <iomanip>      //biblioteka do manipulowanie formatem wyswietlania liczb/napisow
using namespace Eigen;
using namespace std;

//deklaracja sta³ych d³ugoœci cz³onów manipulatora
double d1 = 0.33;   // [m]
double a1 = 0.05;  // [m]
double a2 = 0.33;   // [m]
double a3 = 0.035;  // [m]   
double d2 = 0.335;  // [m]
double d3 = 0.08;   // [m]

//funkcja do tworzenia macierzy transofrmuj¹cej do ka¿dego kolejnego uk³adu T(i->i-1) - zgodnie z notacj¹ Denavita - Hartenberga
Matrix4d DHTransformation(double theta, double d, double a, double phi) {
    Matrix4d transformation;

    transformation << cos(theta), -cos(phi) * sin(theta), sin(phi)* sin(theta), a* cos(theta),
        sin(theta), cos(phi)* cos(theta), -sin(phi) * cos(theta), a* sin(theta),
        0, sin(phi), cos(phi), d,
        0, 0, 0, 1;

    return transformation;
}

//klasa do rozwi¹zywania zadania prostego kinematyki
class ForwardKinematics {

private:

    VectorXd q;                                 //wektor do przechowywania wspó³rzêdnych z³¹czowych, danych w zadaniu prostym
    VectorXd FKresults;                         //wektor do przechowywania wyniku
    double W, P, R, Wdeg, Pdeg, Rdeg;           //zmienne do przechowywania katow po przeliczeniu
    Matrix4d T01, T12, T23, T34, T45, T56, T06; //wszystkie uzywane w klasie macierze 4x4

public:

    //konstruktor klasy FK
    ForwardKinematics(double q1, double q2, double q3, double q4, double q5, double q6)
        :q(6),FKresults(6), W(0.0), P(0.0), R(0.0), Wdeg(0.0), Pdeg(0.0), Rdeg(0.0),
        T01(Matrix4d::Identity()), T12(Matrix4d::Identity()), T23(Matrix4d::Identity()),
        T34(Matrix4d::Identity()), T45(Matrix4d::Identity()), T56(Matrix4d::Identity()),
        T06(Matrix4d::Identity())
    {
        q << q1, q2, q3, q4, q5, q6;
        FKresults.setZero();
    }

    //funkcja sluzaca do rozwiazania FK
    void solveFK() {

        T01 = DHTransformation(q(0), d1, a1, M_PI / 2);         // Macierz T01

        T12 = DHTransformation(q(1), 0.0, a2, 0.0);         // Macierz T12

        T23 = DHTransformation(q(2), 0.0, a3, M_PI / 2);    // Macierz T23

        T34 = DHTransformation(q(3), d2, 0.0, -M_PI / 2);   // Macierz T34

        T45 = DHTransformation(q(4), 0.0, 0.0, M_PI / 2); // Macierz T45

        T56 = DHTransformation(q(5), d3, 0.0, 0.0);          // Macierz T56

        T06 = T01 * T12 * T23 * T34 * T45 * T56;    // wynikowa macierz T06

        //wyliczenie katow W, P, R z macierzy T06. Warunek zabezpieczaj¹cy przed dzieleniem przez 0 (gimbal lock)
        P = asin(-T06(2, 0));
        if (P != M_PI / 2 and P != -M_PI / 2)
        {
            R = atan2(T06(1, 0)/cos(P), T06(0, 0)/ cos(P));
            W = atan2(T06(2, 1)/ cos(P), T06(2, 2)/ cos(P));
        }
        else
        {
            R = 0;
            W = atan2(-T06(1, 2), T06(1, 1));
        }

        //zamiana z radianow na katy
        Wdeg = W * 180 / M_PI;
        Pdeg = P * 180 / M_PI;
        Rdeg = R * 180 / M_PI;
    }

    //funkcja drukujaca macierz T06
    void printT06matrix() {
        cout << "Transformation matrix T06:" << endl << T06 << endl;
    }

    //funkcja drukujaca wyniki FK - polozenia x, y, z [m] i katy W, P, R [deg]
    void printFKResults() {
        cout << endl << "Forward kinematics results:";
        cout << endl << "x = " << T06(0, 3);
        cout << endl << "y = " << T06(1, 3);
        cout << endl << "z = " << T06(2, 3);
        cout << endl;
        cout << endl << "W = " << Wdeg;
        cout << endl << "P = " << Pdeg;
        cout << endl << "R = " << Rdeg;
        cout << endl;
    }

    //funkcja zmieniaca jednostke danych wejsciowych ze stopni na radiany
    void qDegToRad() {
        q = q * (M_PI / 180);
    }

    //funkcja zwracajaca wynik FK w postaci wektora [x y z W P R] (m, deg)
    VectorXd saveFKResults() {
        FKresults << T06(0, 3), T06(1, 3), T06(2, 3), Wdeg, Pdeg, Rdeg;
        return FKresults;
    }

};

//klasa do rozwi¹zywania zadania odwrotnego kinematyki
class InverseKinematics {

private:

    VectorXd position;     //wektor do przechowywania po³o¿enia koñcówki robota w postaci [x y z W P R], dany w zadaniu odwrotnym
    double W, P, R, xc, yc, zc;  //zmienne do przechowywania pozycji
    double beta1front,beta1back, beta2front,beta2back, b, cfront,cback, psifront,psiback, gamma1;  //zmienne do przechowywania wartosci posrednich potrzebnych przy rozwiazywaniu zadania odwrotnego kinematyki
    Matrix4d T06, T01frontdownnoflip, T12frontdownnoflip, T23frontdownnoflip, T03frontdownnoflip, T01frontupnoflip, T12frontupnoflip, T23frontupnoflip, T03frontupnoflip; //wszystkie uzywane w klasie macierze 4x4
    Matrix4d T01backdownnoflip, T12backdownnoflip, T23backdownnoflip, T03backdownnoflip, T01backupnoflip, T12backupnoflip, T23backupnoflip, T03backupnoflip;
    Matrix3d R06, R03frontdownnoflip, R03invfrontdownnoflip, R36frontdownnoflip, R03frontupnoflip, R03invfrontupnoflip, R36frontupnoflip;  //wszystkie uzywane w klasie macierze 3x3
    Matrix3d R03backdownnoflip, R03invbackdownnoflip, R36backdownnoflip, R03backupnoflip, R03invbackupnoflip, R36backupnoflip;

    VectorXd IKresultsfrontdownnoflip; //wektory przechowuj¹ce wynik zadania odwrotnego kinematyki w kazdym mozliwym u³ozeniu w postaci: [theta1 theta2 theta3 theta4 theta5 theta6] (rad)
    VectorXd IKresultsfrontupnoflip; 
    VectorXd IKresultsbackdownnoflip;  
    VectorXd IKresultsbackupnoflip;
    VectorXd IKresultsfrontdownflip;
    VectorXd IKresultsfrontupflip;
    VectorXd IKresultsbackdownflip;
    VectorXd IKresultsbackupflip;

    MatrixXd IKresults; //macierz przechowujaca wszystkie wyniki

public:

    //konstruktor klasy IK
    InverseKinematics(double x, double y, double z, double Wdeg, double Pdeg, double Rdeg)
        :position(6), IKresultsfrontdownnoflip(6), IKresultsfrontupnoflip(6), IKresultsbackdownnoflip(6), IKresultsbackupnoflip(6), IKresultsfrontdownflip(6), IKresultsfrontupflip(6), IKresultsbackdownflip(6), IKresultsbackupflip(6), IKresults(8,6),
        W(0.0), P(0.0), R(0.0), xc(0.0), yc(0.0), zc(0.0), beta1front(0.0), beta1back(0.0), beta2front(0.0), beta2back(0.0), b(0.0), cfront(0.0),cback(0.0), psifront(0.0), gamma1(0.0), psiback(0.0),
        T06(Matrix4d::Identity()), T01frontdownnoflip(Matrix4d::Identity()), T01backdownnoflip(Matrix4d::Identity()), T12frontdownnoflip(Matrix4d::Identity()), T12backdownnoflip(Matrix4d::Identity()),
        T23frontdownnoflip(Matrix4d::Identity()), T23backdownnoflip(Matrix4d::Identity()), T03frontdownnoflip(Matrix4d::Identity()), T03backdownnoflip(Matrix4d::Identity()), R06(Matrix3d::Identity()),
        R03frontdownnoflip(Matrix3d::Identity()), R03backdownnoflip(Matrix3d::Identity()), R03invfrontdownnoflip(Matrix3d::Identity()), R03invbackdownnoflip(Matrix3d::Identity()), R36frontdownnoflip(Matrix3d::Identity()), R36backdownnoflip(Matrix3d::Identity()),
        T01frontupnoflip(Matrix4d::Identity()), T01backupnoflip(Matrix4d::Identity()), T12frontupnoflip(Matrix4d::Identity()), T12backupnoflip(Matrix4d::Identity()), T23frontupnoflip(Matrix4d::Identity()), T23backupnoflip(Matrix4d::Identity()),
        T03frontupnoflip(Matrix4d::Identity()), T03backupnoflip(Matrix4d::Identity()), R03frontupnoflip(Matrix3d::Identity()), R03backupnoflip(Matrix3d::Identity()), R03invfrontupnoflip(Matrix3d::Identity()), R03invbackupnoflip(Matrix3d::Identity()),
        R36frontupnoflip(Matrix3d::Identity()), R36backupnoflip(Matrix3d::Identity())
    {
        position << x, y, z, Wdeg, Pdeg, Rdeg;
        IKresultsfrontdownnoflip.setZero();
        IKresultsfrontupnoflip.setZero();
        IKresultsbackdownnoflip.setZero();
        IKresultsbackupnoflip.setZero();
        IKresultsfrontdownflip.setZero();
        IKresultsfrontupflip.setZero();
        IKresultsbackdownflip.setZero();
        IKresultsbackupflip.setZero();
        W = Wdeg;
        P = Pdeg;
        R = Rdeg;
    }


    //funkcja sluzaca do rozwiazania IK
    void solveIK() {

        //macierz T06 stworzona na podstawie danych do zadania odwrotnego
        T06 << cos(R) * cos(P), -sin(R) * cos(W) + cos(R) * sin(P) * sin(W), sin(R)* sin(W) + cos(R) * sin(P) * cos(W), position(0),
            sin(R)* cos(P), cos(R)* cos(W) + sin(R) * sin(P) * sin(W), -cos(R) * sin(W) + sin(R) * sin(P) * cos(W), position(1),
            -sin(P), cos(P)* sin(W), cos(P)* cos(W), position(2),
            0, 0, 0, 1;

        //pozycja x,y,z kiœci sferycznej (wrist centre)
        xc = T06(0, 3) - T06(0, 2) * d3;
        yc = T06(1, 3) - T06(1, 2) * d3;
        zc = T06(2, 3) - T06(2, 2) * d3;

        //obliczenia dla k¹tów wynikaj¹cych z geometrii: theta1, theta2, theta3

        //theta1 w kazdym mozliwym u³ozeniu
        IKresultsfrontdownnoflip(0) = atan2(yc, xc); 
        IKresultsfrontupnoflip(0) = IKresultsfrontdownnoflip(0);
        IKresultsfrontdownflip(0) = IKresultsfrontdownnoflip(0);
        IKresultsfrontupflip(0) = IKresultsfrontdownnoflip(0);
        if (atan2(yc, xc) < 0)  //warunek wynikajacy z mozliwosc ruchu pierwszego silnika (+-180deg)
        {
            IKresultsbackdownnoflip(0) = atan2(yc, xc) + M_PI; 
            IKresultsbackupnoflip(0) = IKresultsbackdownnoflip(0);
            IKresultsbackdownflip(0) = IKresultsbackdownnoflip(0);
            IKresultsbackupflip(0) = IKresultsbackdownnoflip(0);
        }
        else 
        {
            IKresultsbackdownnoflip(0) = atan2(yc, xc) - M_PI;
            IKresultsbackupnoflip(0) = IKresultsbackdownnoflip(0);
            IKresultsbackdownflip(0) = IKresultsbackdownnoflip(0);
            IKresultsbackupflip(0) = IKresultsbackdownnoflip(0);
        }
        
        //obliczenia opisane w raporcie
        if ((sqrt(xc * xc + yc * yc) - a1) > 0)
        {
            beta1front = atan2(fabs(zc - d1), fabs(sqrt(xc * xc + yc * yc) - a1));
            beta1back = M_PI - atan2(fabs(zc - d1), sqrt(xc * xc + yc * yc) + a1);
        }
        else
        {
            beta1front = M_PI - atan2(fabs(zc - d1), fabs(sqrt(xc * xc + yc * yc) - a1));
            beta1back = M_PI - atan2(fabs(zc - d1), sqrt(xc * xc + yc * yc) + a1);
        }
        b = sqrt(d2 * d2 + a3 * a3);
        cfront = sqrt((sqrt(xc * xc + yc * yc) - a1) * (sqrt(xc * xc + yc * yc) - a1) + (zc - d1) * (zc - d1));
        cback = sqrt((sqrt(xc * xc + yc * yc) + a1) * (sqrt(xc * xc + yc * yc) + a1) + (zc - d1) * (zc - d1));
        beta2front = acos((a2 * a2 + cfront * cfront - b * b) / (2 * a2 * cfront));
        beta2back = acos((a2 * a2 + cback * cback - b * b) / (2 * a2 * cback));

        //theta2 w kazdym mozliwym u³ozeniu
        if (zc > d1)
        {
            IKresultsfrontdownnoflip(1) = beta1front - beta2front;
            IKresultsfrontupnoflip(1) = beta1front + beta2front;
            IKresultsfrontdownflip(1) = beta1front - beta2front;
            IKresultsfrontupflip(1) = beta1front + beta2front;
            IKresultsbackdownnoflip(1) = beta1back - beta2back; 
            IKresultsbackupnoflip(1) = beta1back + beta2back;
            IKresultsbackdownflip(1) = beta1back - beta2back; 
            IKresultsbackupflip(1) = beta1back + beta2back;
        }
        else
        {
            IKresultsfrontdownnoflip(1) = -beta1front - beta2front;
            IKresultsfrontupnoflip(1) = -beta1front + beta2front;
            IKresultsfrontdownflip(1) = -beta1front - beta2front;
            IKresultsfrontupflip(1) = -beta1front + beta2front;
            IKresultsbackdownnoflip(1) = -beta1back - beta2back;
            IKresultsbackupnoflip(1) = -beta1back + beta2back;
            IKresultsbackdownflip(1) = -beta1back - beta2back;
            IKresultsbackupflip(1) = -beta1back + beta2back;
        }

        //warunki if wynikajace z mozliwosci ruchu drugiego silnika
        if (IKresultsbackdownnoflip(1) < -0.96)
        {
            IKresultsbackdownnoflip(1) = IKresultsbackdownnoflip(1) + 2 * M_PI;
            IKresultsbackdownflip(1) = IKresultsbackdownnoflip(1);
        }
        if (IKresultsbackupnoflip(1) < -0.96)
        {
            IKresultsbackupnoflip(1) = IKresultsbackupnoflip(1) + 2 * M_PI;
            IKresultsbackupflip(1) = IKresultsbackupnoflip(1);
        }

        //obliczenia opisane w raporcie
        gamma1 = atan2(a3, d2);
        psifront = acos((a2 * a2 + b * b - cfront * cfront) / (2 * a2 * b));
        psiback = acos((a2 * a2 + b * b - cback * cback) / (2 * a2 * b));

        //theta3 w kazdym mozliwym u³ozeniu
        IKresultsfrontdownnoflip(2) = 3*M_PI/2 - psifront - gamma1; 
        IKresultsfrontupnoflip(2) = -M_PI/2 + psifront - gamma1;
        IKresultsfrontdownflip(2) = 3 * M_PI / 2 - psifront - gamma1;
        IKresultsfrontupflip(2) = -M_PI / 2 + psifront - gamma1;
        IKresultsbackdownnoflip(2) = 3 * M_PI / 2 - psiback - gamma1;
        IKresultsbackupnoflip(2) = -M_PI / 2 + psiback - gamma1;
        IKresultsbackdownflip(2) = 3 * M_PI / 2 - psiback - gamma1; 
        IKresultsbackupflip(2) = -M_PI / 2 + psiback - gamma1;

        //pozosta³e 3 k¹ty: R36 = R03^-1 * R06
        R06 = T06.block<3, 3>(0, 0); //wyodrebnienie z macierzy T06 czêsci zawieraj¹cej R06

        // Macierz T01 w kazdym mozliwym u³ozeniu
        T01frontdownnoflip = DHTransformation(IKresultsfrontdownnoflip(0), d1, a1, M_PI / 2);     
        T01frontupnoflip = DHTransformation(IKresultsfrontupnoflip(0), d1, a1, M_PI / 2);
        T01backdownnoflip = DHTransformation(IKresultsbackdownnoflip(0), d1, a1, M_PI / 2);     
        T01backupnoflip = DHTransformation(IKresultsbackupnoflip(0), d1, a1, M_PI / 2);
        // Macierz T12  w kazdym mozliwym u³ozeniu
        T12frontdownnoflip = DHTransformation(IKresultsfrontdownnoflip(1), 0.0, a2, 0.0);       //tworzymy baze do stworzenia macierzy T03, a pozniej wyodrebnienia z niej R03
        T12frontupnoflip = DHTransformation(IKresultsfrontupnoflip(1), 0.0, a2, 0.0);
        T12backdownnoflip = DHTransformation(IKresultsbackdownnoflip(1), 0.0, a2, 0.0);     
        T12backupnoflip = DHTransformation(IKresultsbackupnoflip(1), 0.0, a2, 0.0);
        // Macierz T23 w kazdym mozliwym u³ozeniu
        T23frontdownnoflip = DHTransformation(IKresultsfrontdownnoflip(2), 0.0, a3, M_PI / 2); 
        T23frontupnoflip = DHTransformation(IKresultsfrontupnoflip(2), 0.0, a3, M_PI / 2);
        T23backdownnoflip = DHTransformation(IKresultsbackdownnoflip(2), 0.0, a3, M_PI / 2);
        T23backupnoflip = DHTransformation(IKresultsbackupnoflip(2), 0.0, a3, M_PI / 2);

        //macierz T03 w kazdym mozliwym u³ozeniu
        T03frontdownnoflip = T01frontdownnoflip * T12frontdownnoflip * T23frontdownnoflip; 
        T03frontupnoflip = T01frontupnoflip * T12frontupnoflip * T23frontupnoflip;
        T03backdownnoflip = T01backdownnoflip * T12backdownnoflip * T23backdownnoflip;
        T03backupnoflip = T01backupnoflip * T12backupnoflip * T23backupnoflip;
        //wyodrebnienie z macierzy T03 czêsci zawieraj¹cej R03
        R03frontdownnoflip = T03frontdownnoflip.block<3, 3>(0, 0); 
        R03frontupnoflip = T03frontupnoflip.block<3, 3>(0, 0);
        R03backdownnoflip = T03backdownnoflip.block<3, 3>(0, 0);
        R03backupnoflip = T03backupnoflip.block<3, 3>(0, 0);
        //odwrocenie macierzy R03
        R03invfrontdownnoflip = R03frontdownnoflip.inverse(); 
        R03invfrontupnoflip = R03frontupnoflip.inverse();
        R03invbackdownnoflip = R03backdownnoflip.inverse();
        R03invbackupnoflip = R03backupnoflip.inverse();

        //macierz R36, z ktorej wyciagniemy theta4, theta5 i theta6
        R36frontdownnoflip = R03invfrontdownnoflip * R06;
        R36frontupnoflip = R03invfrontupnoflip * R06;
        R36backdownnoflip = R03invbackdownnoflip * R06;
        R36backupnoflip = R03invbackupnoflip * R06;


        //theta5 w kazdym mozliwym u³ozeniu front elbow-down
        if (R36frontdownnoflip(2, 2) == 1)     //przy wartoœci 1 acos nie dzia³a³ (?), st¹d if
        {
            IKresultsfrontdownnoflip(4) = 0; 
            IKresultsfrontdownflip(4) = 0;
        }
        else
        {
            IKresultsfrontdownnoflip(4) = acos(R36frontdownnoflip(2, 2));
            IKresultsfrontdownflip(4) = -acos(R36frontdownnoflip(2, 2)); 
        }

        //theta4 i theta6 w kazdym mozliwym u³ozeniu front elbow-dwon
        if (IKresultsfrontdownnoflip(4) != 0 and IKresultsfrontdownnoflip(4) != M_PI)
        {
            IKresultsfrontdownnoflip(3) = atan2(R36frontdownnoflip(1, 2)/sin(IKresultsfrontdownnoflip(4)), R36frontdownnoflip(0, 2)/ sin(IKresultsfrontdownnoflip(4)));
            IKresultsfrontdownnoflip(5) = atan2(R36frontdownnoflip(2, 1)/ sin(IKresultsfrontdownnoflip(4)), -R36frontdownnoflip(2, 0)/ sin(IKresultsfrontdownnoflip(4)));
            IKresultsfrontdownflip(3) = atan2(R36frontdownnoflip(1, 2) / sin(IKresultsfrontdownnoflip(4)), R36frontdownnoflip(0, 2) / sin(IKresultsfrontdownnoflip(4))) - M_PI;
            IKresultsfrontdownflip(5) = atan2(R36frontdownnoflip(2, 1) / sin(IKresultsfrontdownnoflip(4)), -R36frontdownnoflip(2, 0) / sin(IKresultsfrontdownnoflip(4))) + M_PI;
        }
        else
        {
            IKresultsfrontdownnoflip(3) = 0;
            IKresultsfrontdownnoflip(5) = atan2(R36frontdownnoflip(1, 0), R36frontdownnoflip(1, 1));
            IKresultsfrontdownflip(3) = -M_PI;
            IKresultsfrontdownflip(5) = atan2(R36frontdownnoflip(1, 0), R36frontdownnoflip(1, 1))+M_PI;
        }

       
        //theta5 w kazdym mozliwym u³ozeniu front elbow-up
        if (R36frontupnoflip(2, 2) == 1)        //przy wartoœci 1 acos nie dzia³a³ (?), st¹d if
        {
            IKresultsfrontupnoflip(4) = 0;
            IKresultsfrontupflip(4) = 0;
        }
        else
        {
            IKresultsfrontupnoflip(4) = acos(R36frontupnoflip(2, 2));
            IKresultsfrontupflip(4) = -acos(R36frontupnoflip(2, 2));
        }

        //theta4 i theta6 w kazdym mozliwym u³ozeniu front elbow-up
        if (IKresultsfrontupnoflip(4) != 0 and IKresultsfrontupnoflip(4) != M_PI)
        {
            IKresultsfrontupnoflip(3) = atan2(R36frontupnoflip(1, 2)/ sin(IKresultsfrontupnoflip(4)), R36frontupnoflip(0, 2)/ sin(IKresultsfrontupnoflip(4)));  
            IKresultsfrontupnoflip(5) = atan2(R36frontupnoflip(2, 1)/ sin(IKresultsfrontupnoflip(4)), -R36frontupnoflip(2, 0)/ sin(IKresultsfrontupnoflip(4))); 
            IKresultsfrontupflip(3) = atan2(R36frontupnoflip(1, 2) / sin(IKresultsfrontupnoflip(4)), R36frontupnoflip(0, 2) / sin(IKresultsfrontupnoflip(4)))-M_PI;  
            IKresultsfrontupflip(5) = atan2(R36frontupnoflip(2, 1) / sin(IKresultsfrontupnoflip(4)), -R36frontupnoflip(2, 0) / sin(IKresultsfrontupnoflip(4)))+M_PI;
        }
        else
        {
            IKresultsfrontupnoflip(3) = 0; 
            IKresultsfrontupnoflip(5) = atan2(R36frontupnoflip(1, 0), R36frontupnoflip(1, 1)); 
            IKresultsfrontupflip(3) = -M_PI; 
            IKresultsfrontupflip(5) = atan2(R36frontupnoflip(1, 0), R36frontupnoflip(1, 1))+M_PI; 
        }

        //to samo dla back

        //theta5 w kazdym mozliwym u³ozeniu back elbow-down
        if (R36backdownnoflip(2, 2) == 1)     //przy wartoœci 1 acos nie dzia³a³ (?), st¹d if
        {
            IKresultsbackdownnoflip(4) = 0;
            IKresultsbackdownflip(4) = 0;
        }
        else
        {
            IKresultsbackdownnoflip(4) = acos(R36backdownnoflip(2, 2));
            IKresultsbackdownflip(4) = -acos(R36backdownnoflip(2, 2));
        }

        //theta4 i theta6 w kazdym mozliwym u³ozeniu back elbow-dwon
        if (IKresultsbackdownnoflip(4) != 0 and IKresultsbackdownnoflip(4) != M_PI)
        {
            IKresultsbackdownnoflip(3) = atan2(R36backdownnoflip(1, 2) / sin(IKresultsbackdownnoflip(4)), R36backdownnoflip(0, 2) / sin(IKresultsbackdownnoflip(4)));
            IKresultsbackdownnoflip(5) = atan2(R36backdownnoflip(2, 1) / sin(IKresultsbackdownnoflip(4)), -R36backdownnoflip(2, 0) / sin(IKresultsbackdownnoflip(4)));
            IKresultsbackdownflip(3) = atan2(R36backdownnoflip(1, 2) / sin(IKresultsbackdownnoflip(4)), R36backdownnoflip(0, 2) / sin(IKresultsbackdownnoflip(4)))-M_PI;
            IKresultsbackdownflip(5) = atan2(R36backdownnoflip(2, 1) / sin(IKresultsbackdownnoflip(4)), -R36backdownnoflip(2, 0) / sin(IKresultsbackdownnoflip(4)))+M_PI;
        }
        else
        {
            IKresultsbackdownnoflip(3) = 0;
            IKresultsbackdownnoflip(5) = atan2(R36backdownnoflip(1, 0), R36backdownnoflip(1, 1));
            IKresultsbackdownflip(3) = -M_PI;
            IKresultsbackdownflip(5) = atan2(R36backdownnoflip(1, 0), R36backdownnoflip(1, 1))+M_PI;
        }

        //theta5 w kazdym mozliwym u³ozeniu back elbow-up
        if (R36backupnoflip(2, 2) == 1)        //przy wartoœci 1 acos nie dzia³a³ (?), st¹d if
        {
            IKresultsbackupnoflip(4) = 0;
            IKresultsbackupflip(4) = 0;
        }
        else
        {
            IKresultsbackupnoflip(4) = acos(R36backupnoflip(2, 2));
            IKresultsbackupflip(4) = -acos(R36backupnoflip(2, 2)); 
        }

        //theta4 i theta6 w kazdym mozliwym u³ozeniu back elbow-up
        if (IKresultsbackupnoflip(4) != 0 and IKresultsbackupnoflip(4) != M_PI)
        {
            IKresultsbackupnoflip(3) = atan2(R36backupnoflip(1, 2) / sin(IKresultsbackupnoflip(4)), R36backupnoflip(0, 2) / sin(IKresultsbackupnoflip(4)));  
            IKresultsbackupnoflip(5) = atan2(R36backupnoflip(2, 1) / sin(IKresultsbackupnoflip(4)), -R36backupnoflip(2, 0) / sin(IKresultsbackupnoflip(4))); 
            IKresultsbackupflip(3) = atan2(R36backupnoflip(1, 2) / sin(IKresultsbackupnoflip(4)), R36backupnoflip(0, 2) / sin(IKresultsbackupnoflip(4)))-M_PI;  
            IKresultsbackupflip(5) = atan2(R36backupnoflip(2, 1) / sin(IKresultsbackupnoflip(4)), -R36backupnoflip(2, 0) / sin(IKresultsbackupnoflip(4)))+M_PI; 
        }
        else
        {
            IKresultsbackupnoflip(3) = 0; 
            IKresultsbackupnoflip(5) = atan2(R36backupnoflip(1, 0), R36backupnoflip(1, 1));  
            IKresultsbackupflip(3) = -M_PI; 
            IKresultsbackupflip(5) = atan2(R36backupnoflip(1, 0), R36backupnoflip(1, 1))+M_PI;  
        }
    }

    //funkcja drukujaca wyniki IK -  wspolrzedne zlaczowe q1, q2, q3, q4, q5, q6
    void printIKResults() {
        cout << endl << setw(57) << left << "Inverse kinematics results front elbow-down no-flip:" << "Inverse kinematics results back elbow-down no-flip:";
        cout << endl << "theta1 = " << setw(48) << left << IKresultsfrontdownnoflip(0) << "theta1 = " << IKresultsbackdownnoflip(0);
        cout << endl << "theta2 = " << setw(48) << left << IKresultsfrontdownnoflip(1) << "theta2 = " << IKresultsbackdownnoflip(1);
        cout << endl << "theta3 = " << setw(48) << left << IKresultsfrontdownnoflip(2) << "theta3 = " << IKresultsbackdownnoflip(2);
        cout << endl << "theta4 = " << setw(48) << left << IKresultsfrontdownnoflip(3) << "theta4 = " << IKresultsbackdownnoflip(3);
        cout << endl << "theta5 = " << setw(48) << left << IKresultsfrontdownnoflip(4) << "theta5 = " << IKresultsbackdownnoflip(4);
        cout << endl << "theta6 = " << setw(48) << left << IKresultsfrontdownnoflip(5) << "theta6 = " << IKresultsbackdownnoflip(5);
        cout << endl << endl << setw(57) << left << "Inverse kinematics results front elbow-down flip:" << "Inverse kinematics results back elbow-down flip:";
        cout << endl << "theta1 = " << setw(48) << left << IKresultsfrontdownflip(0) << "theta1 = " << IKresultsbackdownflip(0);
        cout << endl << "theta2 = " << setw(48) << left << IKresultsfrontdownflip(1) << "theta2 = " << IKresultsbackdownflip(1);
        cout << endl << "theta3 = " << setw(48) << left << IKresultsfrontdownflip(2) << "theta3 = " << IKresultsbackdownflip(2);
        cout << endl << "theta4 = " << setw(48) << left << IKresultsfrontdownflip(3) << "theta4 = " << IKresultsbackdownflip(3);
        cout << endl << "theta5 = " << setw(48) << left << IKresultsfrontdownflip(4) << "theta5 = " << IKresultsbackdownflip(4);
        cout << endl << "theta6 = " << setw(48) << left << IKresultsfrontdownflip(5) << "theta6 = " << IKresultsbackdownflip(5);
        cout << endl << endl << setw(57) << left << "Inverse kinematics results front elbow-up no-flip:" << "Inverse kinematics results back elbow-up no-flip:";
        cout << endl << "theta1 = " << setw(48) << left << IKresultsfrontupnoflip(0) << "theta1 = " << IKresultsbackupnoflip(0);
        cout << endl << "theta2 = " << setw(48) << left << IKresultsfrontupnoflip(1) << "theta2 = " << IKresultsbackupnoflip(1);
        cout << endl << "theta3 = " << setw(48) << left << IKresultsfrontupnoflip(2) << "theta3 = " << IKresultsbackupnoflip(2);
        cout << endl << "theta4 = " << setw(48) << left << IKresultsfrontupnoflip(3) << "theta4 = " << IKresultsbackupnoflip(3);
        cout << endl << "theta5 = " << setw(48) << left << IKresultsfrontupnoflip(4) << "theta5 = " << IKresultsbackupnoflip(4);
        cout << endl << "theta6 = " << setw(48) << left << IKresultsfrontupnoflip(5) << "theta6 = " << IKresultsbackupnoflip(5);
        cout<<endl<<endl<< setw(57) << left << "Inverse kinematics results front elbow-up flip:" << "Inverse kinematics results back elbow-up flip:";
        cout << endl << "theta1 = " << setw(48) << left << IKresultsfrontupflip(0) << "theta1 = " << IKresultsbackupflip(0);
        cout << endl << "theta2 = " << setw(48) << left << IKresultsfrontupflip(1) << "theta2 = " << IKresultsbackupflip(1);
        cout << endl << "theta3 = " << setw(48) << left << IKresultsfrontupflip(2) << "theta3 = " << IKresultsbackupflip(2);
        cout << endl << "theta4 = " << setw(48) << left << IKresultsfrontupflip(3) << "theta4 = " << IKresultsbackupflip(3);
        cout << endl << "theta5 = " << setw(48) << left << IKresultsfrontupflip(4) << "theta5 = " << IKresultsbackupflip(4);
        cout << endl << "theta6 = " << setw(48) << left << IKresultsfrontupflip(5) << "theta6 = " << IKresultsbackupflip(5);
    }

    //funkcja drukujaca wszystkie wartosci posrednie (wykorzystywana do analizy procedury i szukania bledow)
    void printIntermediateValues(){
        cout << endl;
        cout << endl << "xc = " << xc;
        cout << endl << "yc = " << yc;
        cout << endl << "zc = " << zc;
        cout << endl;
        cout << endl << "zc-d1 = " << zc - d1;
        cout << endl << "sqrt(xc*xc + yc*yc) - a1 = " << sqrt(xc * xc + yc * yc) - a1;
        cout << endl << "beta1 = " << beta1back;
        cout << endl;
        cout << endl << "b = " << b;
        cout << endl << "c = " << cfront;
        cout << endl << "(a2 * a2 + c * c - b * b) / (2 * a2 * c) = " << (a2 * a2 + cfront * cfront - b * b) / (2 * a2 * cfront);       
        cout << endl << "beta2 = " << beta2back;
        cout << endl;
        cout << endl << "gamma1 = " << gamma1 * 180 / M_PI;
        cout << endl;   
        cout << endl << "psi = " << psifront * 180 / M_PI;
        cout << endl;
    }

    //funkcja zmieniaca jednostke danych wejsciowych ze stopni na radiany
    void degToRad() {
        W = position(3) * M_PI / 180;
        P = position(4) * M_PI / 180;
        R = position(5) * M_PI / 180;
    }

    //funkcja zwracajaca wynik IK w postaci macierzy 8x6 z 8 mozliwymi rozwiazaniami 8x[th1 th2 th3 th4 th5 th6] (rad)
    MatrixXd saveIKResults() {
        IKresults << IKresultsfrontdownnoflip(0), IKresultsfrontdownnoflip(1), IKresultsfrontdownnoflip(2), IKresultsfrontdownnoflip(3), IKresultsfrontdownnoflip(4), IKresultsfrontdownnoflip(5),
            IKresultsbackdownnoflip(0), IKresultsbackdownnoflip(1), IKresultsbackdownnoflip(2), IKresultsbackdownnoflip(3), IKresultsbackdownnoflip(4), IKresultsbackdownnoflip(5),
            IKresultsfrontdownflip(0), IKresultsfrontdownflip(1), IKresultsfrontdownflip(2), IKresultsfrontdownflip(3), IKresultsfrontdownflip(4), IKresultsfrontdownflip(5),
            IKresultsbackdownflip(0), IKresultsbackdownflip(1), IKresultsbackdownflip(2), IKresultsbackdownflip(3), IKresultsbackdownflip(4), IKresultsbackdownflip(5),
            IKresultsfrontupnoflip(0), IKresultsfrontupnoflip(1), IKresultsfrontupnoflip(2), IKresultsfrontupnoflip(3), IKresultsfrontupnoflip(4), IKresultsfrontupnoflip(5),
            IKresultsbackupnoflip(0), IKresultsbackupnoflip(1), IKresultsbackupnoflip(2), IKresultsbackupnoflip(3), IKresultsbackupnoflip(4), IKresultsbackupnoflip(5),
            IKresultsfrontupflip(0), IKresultsfrontupflip(1), IKresultsfrontupflip(2), IKresultsfrontupflip(3), IKresultsfrontupflip(4), IKresultsfrontupflip(5),
            IKresultsbackupflip(0), IKresultsbackupflip(1), IKresultsbackupflip(2), IKresultsbackupflip(3), IKresultsbackupflip(4), IKresultsbackupflip(5);
        return IKresults;
    }

};

//funkcja wywo³uj¹ca wszystkie procedury dla zadania prostego
void FK(VectorXd& FKin, VectorXd& FKout)
{
    ForwardKinematics fk1(FKin(0), FKin(1), FKin(2), FKin(3), FKin(4), FKin(5));        //tworzenie obiektu klasy FK

    fk1.qDegToRad();                 //uzyj jesli dane wejsciowe w s¹ w stopniach
    fk1.solveFK();                  //uzyj zeby rozwiazac kinematyke prosta
    fk1.printT06matrix();           //uzyj zeby wydrukowac macierz T06
    fk1.printFKResults();           //uzyj zeby wydrukowac wyniki zadania prostego
    FKout = fk1.saveFKResults();    //uzyj zeby zapisac wyniki do wektora na zewnatrz klasy
}

//funkcja wywo³uj¹ca wszystkie procedury dla zadania odwrotnego
void IK(VectorXd& IKin, MatrixXd& IKout)
{
    InverseKinematics ik1(IKin(0), IKin(1), IKin(2), IKin(3), IKin(4), IKin(5));        //tworzenie obiektu klasy IK

    ik1.degToRad();                     //uzyj jesli dane wejsciowe s¹ w stopniach
    ik1.solveIK();                      //uzyj zeby rozwiazac kinematyke odwrotna
    ik1.printIKResults();             //uzyj zeby wydrukowac wyniki zadania odwrotnego
    //ik1.printIntermediateValues();    //uzyj zeby wydrukowac wyniki posrednie (dla analizy, sprawdzenia)
    IKout = ik1.saveIKResults();        //uzyj zeby zapisac wyniki do macierzy na zewnatrz klasy
}

//funkcja do uruchomienia "interfejsu" w celu testowania programu
void test()
{
    int choice;         //zmienna przechowuj¹ca odpowiedz uzytkownika

    //wektory i macierz dzia³ajacê w funkcji test()
    VectorXd FKi(6);
    VectorXd IKi(6);
    VectorXd FKo(6);
    MatrixXd IKo(8,6);

    while(1)
    {
        cout << "Kinematics to solve:" << endl << "1 - forward" << endl << "2 - inverse" << endl << "3 - clear the screen"<<endl;
        cin >> choice;
        if (choice == 1)
        {
            cout << endl << "Enter your input set one by one (in degrees):" << endl << "theta1 = ";
            cin >> FKi(0);
            cout << "theta2 = ";
            cin >> FKi(1);
            cout << "theta3 = ";
            cin >> FKi(2);
            cout << "theta4 = ";
            cin >> FKi(3);
            cout << "theta5 = ";
            cin >> FKi(4);
            cout << "theta6 = ";
            cin >> FKi(5);
            cout << endl;
            FK(FKi, FKo);
            cout << endl;
        }
        else if (choice == 2)
        {
            cout << endl << "Enter you input set one by one (in meters and degrees):" << endl << "x = ";
            cin >> IKi(0);
            cout << "y = ";
            cin >> IKi(1);
            cout << "z = ";
            cin >> IKi(2);
            cout << "W = ";
            cin >> IKi(3);
            cout << "P = ";
            cin >> IKi(4);
            cout << "R = ";
            cin >> IKi(5);
            cout << endl;
            IK(IKi, IKo);
            cout << endl<<endl;      
        }
        else if(choice==3)
        {
            system("cls");
        }
        else
        {
            cout << "error" << endl;
        }
    }
}

//funkcja do czytania macierzy z pliku csv z losowo wygenerowanymi konfiguracjami robota
MatrixXd readMatrix(const string& filename, int rows, int cols) {
    MatrixXd matrix(rows, cols);

    //otworzenie pliku
    ifstream file(filename);    
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return matrix;
    }

    //pobieranie danych z pliku i zapisywanie ich do utworzonej macierzy zamieniaj¹c ze string na double 
    string line;
    int row = 0;
    while (getline(file, line) && row < rows) {
        stringstream ss(line);
        string value;
        int col = 0;

        while (getline(ss, value, ',') && col < cols) {
            matrix(row, col) = stod(value);
            col++;
        }
        row++;
    }
    
    file.close(); //zamkniecie pliku
    return matrix;
}

//funkcja przeprowadzajaca test dzialania programu dla pliku z wygenerowanymi konfiguracjami
void testGenereatedConfigurations(VectorXd& FKin, VectorXd& FKout, VectorXd& IKin, MatrixXd& IKout)
{
    string filename = "generatedConfigurations.csv";    //nazwa pliku na którym dzia³amy

    //wymiary zawartej w pliku macierzy danych
    int rows = 10000;
    int cols = 6;

    MatrixXd generatedConfigs = readMatrix(filename, rows, cols);  //pobranie macierzy z pliku z wykorzystaniem stworzonej wczesniej funkcji

    //potrzebne do przeprowadzenia analizy wektory
    VectorXd vector(6);
    VectorXd onesetnorms(8);
    VectorXd allsetsnorms(10000);

    //procedura do przeprowadzenia analizy dzialania programu szerzej opisana w raporcie
    for (int i = 0; i < 10000; i++)
    {
        FKin = generatedConfigs.row(i);     //pobieranie kolejnych wierszy macierzy na wejscie do zadania prostego
        FK(FKin, FKout);    //zad. proste
        IKin = FKout;           //przypisanie wyjscia zad. prostego na wejscie do zadania odwrotnego
        IK(IKin, IKout);    //zad. odwrotne

        //for wynikajacy z faktu ze nalezy sprawdzic, ktore z 8 mozliwych rozwiazan jest aktualnie poszukiwanym
        for (int j = 0; j < 8; j++)
        {
            vector << FKin(0) - IKout(j, 0), FKin(1) - IKout(j, 1), FKin(2) - IKout(j, 2), FKin(3) - IKout(j, 3), FKin(4) - IKout(j, 4), FKin(5) - IKout(j, 5);

            //for uwzgledniajacy przesuniecie rozwiazan o 360deg (flip / no-flip)
            for (int k = 3; k < 6; k++)
            {
                if (vector(k) > 2 * M_PI - 0.01 && vector(k) < 6.3 + 0.01)
                {
                    vector(k) = vector(k) - 2 * M_PI;
                }
                else if (vector(k) < -2 * M_PI + 0.01 && vector(k) > -2 * M_PI - 0.01)
                {
                    vector(k) = vector(k) + 2 * M_PI;
                }
            }
            onesetnorms(j) = vector.norm(); //wektor z normami wszystkich 8 przypadkow
        }
        allsetsnorms(i) = onesetnorms.minCoeff();   //pobranie do ogolnego wektora 1x10000 tylko najmniejszej normy (rozwiazania, ktorego poszukiwano)
    }
    cout << "max = " << allsetsnorms.maxCoeff() << endl;    //maksymalna norma wsrod wszystkich 10000 sprawdzanych konfiguracji
    cout << "mean = " << allsetsnorms.mean() << endl;       //srednia z wszystkich policzonych norm
}


int main()
{
    //inicjalizacja wektorow oraz macierzy do przechowywania danych wejœciowych oraz wyników
    VectorXd FKin(6);
    VectorXd IKin(6);
    VectorXd FKout(6);
    MatrixXd IKout(8,6);

    FKin << 0, 0, 0, 0, 0, 0;         //przykladowy zestaw danych [th1 th2 th3 th4 th5 th6] (deg)
 
    test();             //uzyj, jesli chcesz przetestowac dzia³anie funkcji wprowadzajac dane na bie¿¹co

    FK(FKin, FKout);    //uzyj, zeby rozwiazac zadanie proste od poczatku do konca
    
    IKin = FKout;     //przykladowy zestaw danych [x y z W P R] (m, deg)

    IK(IKin, IKout);    //uzyj, zeby rozwiazac zadanie odwrotne od poczatku do konca

    //testGenereatedConfigurations(FKin, FKout, IKin, IKout);  //uzyj, zeby przeprowadzic test dzialania programu dla calego zestawu wczesniej przygotowanych konfiguracji
                                                                        // w FK nalezy skomentowac zamiane katow na radiany (dane w rad). zaleca sie wylaczenie wyswietlania wynikow w FK i IK
    return 0;
}