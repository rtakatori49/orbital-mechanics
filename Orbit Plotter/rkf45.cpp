// Runge-Kutta-Fehlberg Method
// Ryo Takatori
// 09/23/2020

int main(){
    double h;
    
    k1 = h*f(t,y);
    k2 = h*f(t+(1/4)*h,y+(1/4)*k1);
    k3 = h*f(t+(3/8)*h,y+(3/32)*k1+(9/32)*k2);
    k4 = h*f(t+(12/13)*h,y)
}