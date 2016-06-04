disp('--- (©) Alexei Dergunov, 2009 ---');

//----------------------------------------
disp('----------- Задание 3 -----------');
//----------------------------------------

funcprot(0);
n = 400; disp(n, 'объем выборки');
ax = 2.75; disp(ax, 'мат. ожидание X');
dx = 6.5; disp(dx, 'дисперсия X');
ay = -5.3; disp(ay, 'мат. ожидание Y');
dy = 1.8; disp(dy, 'дисперсия Y');
rxy = 0.11; disp(rxy, 'коэффициент корреляции');
alpha = 0.01; disp(alpha, 'уровень значимости');

function[y] = f(x,a,d) // нормальный закон распределения N(a,d)
  y = (1/sqrt(2*%pi*d)) * exp(-(x-a)^2 / (2*d));
endfunction

// выборки случайных величин X~N(ax,dx) и Y~N(ay,dy)
Z1 = grand(1,n,'nor',0,1);
Z2 = grand(1,n,'nor',0,1);
B = [sqrt(dx),     0               ; ..
     rxy*sqrt(dy), sqrt(dy*(1-rxy))  ..
    ]; // матрица преобразования
X = B(1,1)*Z1 + ax;
Y = B(2,1)*Z1 + B(2,2)*Z2 + ay;


//------------------------------------------------------
disp('');
disp('Задание 3.1');

xsr = mean(X);
disp(xsr, 'выборочное среднее случайной величины X');
disp(abs(ax-xsr),'|MX - x(среднее)|');
ysr = mean(Y);
disp(ysr, 'выборочное среднее случайной величины Y');
disp(abs(ay-ysr),'|MY - y(среднее)|');
sx2 = mean(X^2) - xsr^2;
disp(sx2, 'выборочная дисперсия случайной величины X');
disp(abs(dx-sx2),'|DX - sx^2|');
sy2 = mean(Y^2) - ysr^2;
disp(sy2, 'выборочная дисперсия случайной величины Y');
disp(abs(dy-sy2),'|DY - sy^2|');
xysr=0; for i=1:n, xysr=xysr+X(i)*Y(i); end;
r_xy = (xysr/n - xsr*ysr) / sqrt(sx2*sy2);
disp(r_xy, 'выборочный коэффициент корреляции');
disp(abs(rxy-r_xy),'|rxy - r^xy|');


//------------------------------------------------------
disp('');
disp('Задание 3.2');

k = floor(1 + 3.32*log10(n)) + 1;
disp(k, 'число интервалов');
ux = zeros(1,k+1); // границы интервалов X
uy = zeros(1,k+1); // границы интервалов Y
nux = zeros(1,k); // частоты интервалов X
nuy = zeros(1,k); // частоты интервалов Y
nuxy = zeros(k,k); // частоты прямоугольников

// вычисление границ интервалов
ux(1) = min(X);
ux(k+1) = max(X);
uy(1) = min(Y);
uy(k+1) = max(Y);
dux = (ux(k+1) - ux(1)) / k;
duy = (uy(k+1) - uy(1)) / k;
for i = 2 : k,
  ux(i) = ux(i-1) + dux;
  uy(i) = uy(i-1) + duy;
end;

K = k; // число интервалов X
L = k; // число интервалов Y

function[nux,nuy,nuxy] = countFrequencies(ux,uy)
  // вычисление частот
  K = size(ux,2) - 1;
  L = size(uy,2) - 1;
  nux = zeros(1,K); // частоты интервалов X
  nuy = zeros(1,L); // частоты интервалов Y
  nuxy = zeros(K,L); // частоты прямоугольников
  for i = 1 : n,
    jx=0; jy=0;
    for j = K : -1 : 1, // для Х
      if ux(j)<=X(i) & X(i)<=ux(j+1),
        nux(j) = nux(j)+1;
        jx = j;
        break;
      end;
    end;
    for j = L : -1 : 1, // для Y
      if uy(j)<=Y(i) & Y(i)<=uy(j+1),
        nuy(j) = nuy(j)+1;
        jy = j;
        break;
      end;
    end;
    nuxy(jx,jy) = nuxy(jx,jy) + 1;
  end;
endfunction

[nux,nuy,nuxy] = countFrequencies(ux,uy);

disp(ux, 'границы интервалов X');
disp(nux, 'частоты интервалов X');
disp(nux/n, 'отн. частоты интервалов X');
disp(uy, 'границы интервалов Y');
disp(nuy, 'частоты интервалов Y');
disp(nuy/n, 'отн. частоты интервалов Y');
disp(nuxy, 'частоты прямоугольников');
disp(nuxy/n, 'отн. частоты прямоугольников');

function[ans,ii,jj] = minInRect(nuxy)
  // находит числа меньше 5 в матрице частот
  ans = 5; // минимальное число
  ii=1; jj=1; // строка и столбец
  [K,L] = size(nuxy);
  for i = 1 : K,
    for j = 1 : L,
      if nuxy(i,j)<ans,
        ans = nuxy(i,j);
        ii = i;
        jj = j;
        return;
      end;
    end;
  end;  
endfunction

[min_in_rect,ii,jj] = minInRect(nuxy);
while min_in_rect < 5,
// работает, пока не выполнено условие min >= 5
  if ii==K, ii=K-1; end;
  if jj==L, jj=L-1; end;
  if K > L,
    for i = ii+1 : K,
      ux(i) = ux(i+1);
    end;
    ux = ux(1:K);
    K = K - 1;
  else,
    for j = jj+1 : L,
      uy(j) = uy(j+1);
    end;
    uy = uy(1:L);
    L = L - 1;
  end;
  [nux,nuy,nuxy]=countFrequencies(ux,uy);
  [min_in_rect,ii,jj] = minInRect(nuxy);
end;

disp(K, 'новое число интервалов X');
disp(ux, 'новые границы интервалов X');
disp(nux, 'новые частоты интервалов X');
disp(nux/n, 'новые отн. частоты интервалов X');
disp(L, 'новое число интервалов Y');
disp(uy, 'новые границы интервалов Y');
disp(nuy, 'новые частоты интервалов Y');
disp(nuy/n, 'новые отн. частоты интервалов Y');
disp(nuxy, 'новые частоты прямоугольников');
disp(nuxy/n, 'отн. частоты новых прямоугольников');

px = zeros(1,K); // относительные частоты X
py = zeros(1,L); // относительные частоты Y
pxy = zeros(K,L); // отн. частоты прямоугольников

for i = 1 : K, px(i) = nux(i) / n; end;
for i = 1 : L, py(i) = nuy(i) / n; end;
for i = 1 : K, for j = 1 : L,
  pxy(i,j) = nux(i)*nuy(j) / n^2;
end; end;

// используется критерий проверки значимости rxy
h = 2.576; // порог = (1-alpha)/2 - квантиль
           // распределения Стьюдента S(n-2)
t = (r_xy*sqrt(n-2)) / sqrt(1-r_xy^2); // статистика

// используется критерий независимости hi^2
hh = cdfnor("X", (K-1)*(L-1), 2*(K-1)*(L-1), ..
     1-alpha, alpha); // порог
//hh = cdfnor("X", (K-1)*(L-1), 2*(K-1)*(L-1), ..
//     1-alpha, alpha); // порог
st = 0; 
for i = 1 : K,
  for j = 1 : L,
    st = st + (nuxy(i,j)^2) / (nux(i)*nuy(j));
  end;
end;
st = n*(st - 1); // статистика

disp('используется критерий проверки значимости rxy');
if abs(t) >= h,
  disp([abs(t),h],..
  'гипотеза о независимости отвержена (|t| >= h)');
else
  disp([abs(t),h],..
  'гипотеза о независимости принята (|t| < h)');
end;
disp('используется критерий независимости hi^2');
if abs(st) >= hh,
  disp([abs(st),hh],..
  'гипотеза о независимости отвержена (|t| >= h)');
else
  disp([abs(st),hh],..
  'гипотеза о независимости принята (|t| < h)');
end;


//------------------------------------------------------
disp('');
disp('Задание 3.3');

function [y] = regrYtoX(x) // регрессия Y на X
  y = ysr + r_xy*sqrt(sy2/sx2)*(x-xsr);
endfunction

function [x] = regrXtoY(y) // регрессия X на Y
  x = xsr + r_xy*sqrt(sx2/sy2)*(y-ysr);
endfunction

x = [min(X)-1 : 0.001 : max(X)+1]; // значения X
y = [min(Y)-1 : 0.001 : max(Y)+1]; // значения Y

plot2d(X,Y,0); // значения выборки
plot2d(x,regrYtoX(x),5); // ур-е регрессии Y на X
plot2d(regrXtoY(y),y,5); // ур-е регрессии X на Y

disp('графики построены!');
