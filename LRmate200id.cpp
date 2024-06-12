#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Dense>      //biblioteka umo¿liwiaj¹ca tworzenie oraz wykonywanie dzia³an na macierzach i wektorach
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
        Wdeg = W * (180 / M_PI);
        Pdeg = P * (180 / M_PI);
        Rdeg = R * (180 / M_PI);
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

    //funkcja zwracajaca wynik FK w postaci wektora [x y z W P R] 
    VectorXd saveFKResults() {
        FKresults << T06(0, 3), T06(1, 3), T06(2, 3), Wdeg, Pdeg, Rdeg;
        return FKresults;
    }

};

//klasa do rozwi¹zywania zadania odwrotnego kinematyki
class InverseKinematics {

private:

    VectorXd position;                                                              //wektor do przechowywania po³o¿enia koñcówki robota w postaci [x y z W P R], dany w zadaniu odwrotnym
    double W, P, R, xc, yc, zc;                                                     //zmienne do przechowywania pozycji
    double beta1, beta2, b, c, psi, gamma1, gamma2, gamma;                           //zmienne do przechowywania wartosci posrednich potrzebnych przy rozwiazywaniu zadania odwrotnego kinematyki
    Matrix4d T06, T01down, T12down, T23down, T03down, T01up, T12up, T23up, T03up;   //wszystkie uzywane w klasie macierze 4x4
    Matrix3d R06, R03down, R03invdown, R36down, R03up, R03invup, R36up;             //wszystkie uzywane w klasie macierze 3x3
    VectorXd IKresultsdown;                                                         //wektor przechowuj¹cy wynik zadania odwrotnego kinematyki gdy elbow-down w postaci: [theta1 theta2 theta3 theta4 theta5 theta6]
    VectorXd IKresultsup;                                                           //wektor przechowuj¹cy wynik zadania odwrotnego kinematyki gdy elbow-up w postaci: [theta1 theta2 theta3 theta4 theta5 theta6]

public:

    //konstruktor klasy IK
    InverseKinematics(double x, double y, double z, double Wdeg, double Pdeg, double Rdeg)
        :position(6), IKresultsdown(6), IKresultsup(6), W(0.0), P(0.0), R(0.0), xc(0.0), yc(0.0), zc(0.0),
        beta1(0.0), beta2(0.0), b(0.0), c(0.0), psi(0.0), gamma1(0.0), gamma2(0.0), gamma(0.0),
        T06(Matrix4d::Identity()), T01down(Matrix4d::Identity()), T12down(Matrix4d::Identity()),
        T23down(Matrix4d::Identity()), T03down(Matrix4d::Identity()), R06(Matrix3d::Identity()),
        R03down(Matrix3d::Identity()), R03invdown(Matrix3d::Identity()), R36down(Matrix3d::Identity()), 
        T01up(Matrix4d::Identity()), T12up(Matrix4d::Identity()),  T23up(Matrix4d::Identity()),
        T03up(Matrix4d::Identity()), R03up(Matrix3d::Identity()), R03invup(Matrix3d::Identity()),
        R36up(Matrix3d::Identity())
    {
        position << x, y, z, Wdeg, Pdeg, Rdeg;
        IKresultsdown.setZero();
        IKresultsup.setZero();
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

        //k¹ty wynikaj¹ce z geometrii (IKresults << theta1, theta2, theta3, theta4, theta5, theta6)
        IKresultsdown(0) = atan2(yc, xc); //theta1
        IKresultsup(0) = IKresultsdown(0);
        
        if ((sqrt(xc * xc + yc * yc) - a1) > 0)
        {
            beta1 = atan2(fabs(zc - d1), fabs(sqrt(xc * xc + yc * yc) - a1));
        }
        else
        {
            beta1 = M_PI - atan2(fabs(zc - d1), fabs(sqrt(xc * xc + yc * yc) - a1));
        }
        b = sqrt(d2 * d2 + a3 * a3);
        c = sqrt((sqrt(xc * xc + yc * yc) - a1) * (sqrt(xc * xc + yc * yc) - a1) + (zc - d1) * (zc - d1));
        beta2 = acos((a2 * a2 + c * c - b * b) / (2 * a2 * c));
        
        if (zc > d1)
        {
            IKresultsdown(1) = beta1 - beta2; //theta2
            IKresultsup(1) = beta1 + beta2;
        }
        else
        {
            IKresultsdown(1) = -beta1 - beta2;
            IKresultsup(1) = -beta1 + beta2;
        }

        gamma1 = atan2(a3, d2);
        psi = acos((a2 * a2 + b * b - c * c) / (2 * a2 * b));
        gamma2 = psi - M_PI / 2 - IKresultsdown(1);
        gamma = gamma1 + gamma2;
        IKresultsdown(2) = 3*M_PI/2 - psi - gamma1; //theta3 
        IKresultsup(2) = -M_PI/2 + psi - gamma1; 

        //pozosta³e 3 k¹ty: R36 = R03^-1 * R06
        R06 = T06.block<3, 3>(0, 0); //wyodrebnienie z macierzy T06 czêsci zawieraj¹cej R06

        T01down = DHTransformation(IKresultsdown(0), d1, a1, M_PI / 2);     // Macierz T01
        T01up = DHTransformation(IKresultsup(0), d1, a1, M_PI / 2);
        T12down = DHTransformation(IKresultsdown(1), 0.0, a2, 0.0);      // Macierz T12       tworzymy baze do stworzenia macierzy T03, a pozniej wyodrebnienia z niej R03
        T12up = DHTransformation(IKresultsup(1), 0.0, a2, 0.0);
        T23down = DHTransformation(IKresultsdown(2), 0.0, a3, M_PI / 2); // Macierz T23
        T23up = DHTransformation(IKresultsup(2), 0.0, a3, M_PI / 2);
        T03down = T01down * T12down * T23down; //macierz T03
        T03up = T01up * T12up * T23up;
        R03down = T03down.block<3, 3>(0, 0); //wyodrebnienie z macierzy T03 czêsci zawieraj¹cej R03
        R03up = T03up.block<3, 3>(0, 0);
        R03invdown = R03down.inverse(); //odwrocenie macierzy R03
        R03invup = R03up.inverse();

        R36down = R03invdown * R06; //macierz R36, z ktorej wyciagniemy theta4, theta5 i theta6
        R36up = R03invup * R06;


       
        if (R36down(2, 2) == 1)     //przy wartoœci 1 acos nie dzia³a³ (?), st¹d if
        {
            IKresultsdown(4) = 0;
        }
        else
        {
            IKresultsdown(4) = acos(R36down(2, 2)); //theta5
        }


        if (IKresultsdown(4) != 0 and IKresultsdown(4) != M_PI)
        {
            IKresultsdown(3) = atan2(R36down(1, 2)/sin(IKresultsdown(4)), R36down(0, 2)/ sin(IKresultsdown(4))); //theta4
            IKresultsdown(5) = atan2(R36down(2, 1)/ sin(IKresultsdown(4)), -R36down(2, 0)/ sin(IKresultsdown(4))); //theta6
        }
        else
        {
            IKresultsdown(3) = 0;
            IKresultsdown(5) = atan2(R36down(1, 0), R36down(1, 1));
        }

       
        if (R36up(2, 2) == 1)        //przy wartoœci 1 acos nie dzia³a³ (?), st¹d if
        {
            IKresultsup(4) = 0;
        }
        else
        {
            IKresultsup(4) = acos(R36up(2, 2)); //theta5
        }

        if (IKresultsup(4) != 0 and IKresultsup(4) != M_PI)
        {
            IKresultsup(3) = atan2(R36up(1, 2)/ sin(IKresultsup(4)), R36up(0, 2)/ sin(IKresultsup(4)));  //theta4
            IKresultsup(5) = atan2(R36up(2, 1)/ sin(IKresultsup(4)), -R36up(2, 0)/ sin(IKresultsup(4))); //theta6
        }
        else
        {
            IKresultsup(3) = 0;
            IKresultsup(5) = atan2(R36up(1, 0), R36up(1, 1));
        }

    }

    //funkcja drukujaca wyniki IK -  wspolrzedne zlaczowe q1, q2, q3, q4, q5, q6
    void printIKResults() {
        cout << endl << "Inverse kinematics results option I (elbow-down):";
        cout << endl << "theta1 = " << IKresultsdown(0) * 180 / M_PI;
        cout << endl << "theta2 = " << IKresultsdown(1) * 180 / M_PI;
        cout << endl << "theta3 = " << IKresultsdown(2) * 180 / M_PI;
        cout << endl << "theta4 = " << IKresultsdown(3) * 180 / M_PI << " +-"<<180;
        cout << endl << "theta5 = " << IKresultsdown(4) * 180 / M_PI << " or " << -IKresultsdown(4) * 180 / M_PI;
        cout << endl << "theta6 = " << IKresultsdown(5) * 180 / M_PI << " +-" << 180;

        cout << endl << endl<< "Inverse kinematics results option II (elbow-up):";
        cout << endl << "theta1 = " << IKresultsup(0) * 180 / M_PI;
        cout << endl << "theta2 = " << IKresultsup(1) * 180 / M_PI;
        cout << endl << "theta3 = " << IKresultsup(2) * 180 / M_PI;
        cout << endl << "theta4 = " << IKresultsup(3) * 180 / M_PI << "  +-" << 180;
        cout << endl << "theta5 = " << IKresultsup(4) * 180 / M_PI << "  or " << -IKresultsup(4) * 180 / M_PI;
        cout << endl << "theta6 = " << IKresultsup(5) * 180 / M_PI << "  +-" << 180;
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
        cout << endl << "beta1 = " << beta1 * 180 / M_PI;
        cout << endl;
        cout << endl << "b = " << b;
        cout << endl << "c = " << c;
        cout << endl << "(a2 * a2 + c * c - b * b) / (2 * a2 * c) = " << (a2 * a2 + c * c - b * b) / (2 * a2 * c);       
        cout << endl << "beta2 = " << beta2 * 180 / M_PI;
        cout << endl;
        cout << endl << "gamma1 = " << gamma1 * 180 / M_PI;
        cout << endl;   
        cout << endl << "psi = " << psi * 180 / M_PI;
        cout << endl;
        cout << endl << "gamma = " << gamma * 180 / M_PI;
    }

    //funkcja zmieniaca jednostke danych wejsciowych ze stopni na radiany
    void degToRad() {
        W = position(3) * (M_PI / 180);
        P = position(4) * (M_PI / 180);
        R = position(5) * (M_PI / 180);
    }

    //funkcja zwracajaca wynik IK w postaci wektora [th1 th2 th3 th4 th5 th6] 
    VectorXd saveIKResults() {
        return IKresultsdown * (180 / M_PI);
    }

};

//funkcja wywo³uj¹ca wszystkie procedury dla zadania prostego
void FK(VectorXd FKin, VectorXd FKout)
{
    ForwardKinematics fk1(FKin(0), FKin(1), FKin(2), FKin(3), FKin(4), FKin(5));        //tworzenie obiektu klasy FK

    fk1.qDegToRad();                 //uzyj jesli dane w s¹ w stopniach
    fk1.solveFK();                  //uzyj zeby rozwiazac kinematyke prosta
    //fk1.printT06matrix();           //uzyj zeby wydrukowac macierz T06
    fk1.printFKResults();           //uzyj zeby wydrukowac wyniki zadania prostego
    FKout = fk1.saveFKResults();    //uzyj zeby zapisac wyniki do wektora na zewnatrz klasy
}

//funkcja wywo³uj¹ca wszystkie procedury dla zadania odwrotnego
void IK(VectorXd IKin, VectorXd IKout)
{
    InverseKinematics ik1(IKin(0), IKin(1), IKin(2), IKin(3), IKin(4), IKin(5));        //tworzenie obiektu klasy IK

    ik1.degToRad();                     //uzyj jesli dane w s¹ w stopniach
    ik1.solveIK();                      //uzyj zeby rozwiazac kinematyke odwrotna
    ik1.printIKResults();               //uzyj zeby wydrukowac wyniki zadania odwrotnego
    //ik1.printIntermediateValues();    //uzyj zeby wydrukowac wyniki posrednie (dla analizy, sprawdzenia)
    IKout = ik1.saveIKResults();        //uzyj zeby zapisac wyniki do wektora na zewnatrz klasy
}

//funkcja do uruchomienia "interfejsu" w celu testowania programu
void test()
{
    int choice;         //zmienna przechowuj¹ca odpowiedz uzytkownika

    //wektory dzia³ajacê w funkcji test()
    VectorXd FKi(6);
    VectorXd IKi(6);
    VectorXd FKo(6);
    VectorXd IKo(6);

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

int main()
{
    //inicjalizacja wektorow do przechowywania danych wejœciowych oraz wyników
    VectorXd FKin(6);
    VectorXd IKin(6);
    VectorXd FKout(6);
    VectorXd IKout(6);

    FKin << 32, 100, 100, 23, 65, -58;          //przykladowy zestaw danych [th1 th2 th3 th4 th5 th6] (deg)
    IKin << 0.052, 0.050, 0.2, 43, -65, 71;     //przykladowy zestaw danych [x y z W P R] (m, deg)

    test();             //uzyj, jesli chcesz przetestowac dzia³anie funkcji wprowadzajac dane na bie¿¹co

    FK(FKin, FKout);    //uzyj, zeby rozwiazac zadanie proste od poczatku do konca
    IK(IKin, IKout);    //uzyj, zeby rozwiazac zadanie odwrotne od poczatku do konca

    return 0;
}