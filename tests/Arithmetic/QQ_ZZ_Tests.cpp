/*
 * @file QQ_ZZ_Tests.cpp
 * @author Reginald Lybbert
 * @brief Test program for QQ<ZZ> class for rational numbers and functions
 */


#include <ANTL/Arithmetic/QQ.hpp>


NTL_CLIENT
using namespace ANTL;

int main(){
  
  int failureCount = 0;

  // Base type long

  ZZ a,b,c,d,g;
  
  //get random values for a,b,c,d
  RandomLen(a,128);
  RandomLen(b,128);
  RandomLen(c,128);
  RandomLen(d,128);
  GCD(g,c,d);

  //Test all constructors
  QQ<ZZ> zeroq;
  QQ<ZZ> aq(a);
  QQ<ZZ> bq(b);
  QQ<ZZ> cdq(c,d);
  QQ<ZZ> aCopyq(aq); 
  QQ<ZZ> rq(a,a*b);

  if(IsZero(zeroq.getNumerator()) && IsOne(zeroq.getDenominator())){
      cout << "Empty Constructor: Success" << endl;
  }
  else{
      cout << "Empty Constructor: Failure" << endl;
      failureCount++;
  }

  if(aq.getN() == a && IsOne(aq.getD())){
      cout << "Constructor from one value: Success" << endl; 
  }
  else{
      cout << "Constructor from one value: Failure" << endl;
      failureCount++;
  }

  if(cdq.getN() == c/g && cdq.getD() == d/g){
     cout << "Constructor from two values: Success" << endl;
  }
  else{
     cout << "Constructor from two values: Failure" << endl;
     failureCount++;
  }

  if(aCopyq.getN() == a && IsOne(aq.getD())){
     cout << "Constructor from QQ: Success" << endl;
  }
  else{
     cout << "Constructor from QQ: Failure" << endl;
     failureCount++;
  }

  if(IsOne(rq.getN()) && rq.getD() == b){
      cout << "Constructor with reduction: Success" << endl; 
  }
  else{
      cout << "Constructor with reduction: Failure" << endl;
      failureCount++;
  }

  //test assignment

  cdq.setNumerator(d+1);
  aq.setN(a+1);
  if(cdq.getN() == d+1 && cdq.getD() == d/g){
     cout << "setNumerator: Success" << endl; 
  }
  else{
     cout << "setNumerator: Failure" << endl;
     failureCount++;
  }
 
  if(aq.getN() == a+1 && IsOne(aq.getD())){
     cout << "setN: Success" << endl;
  }
  else{
     cout << "setN: Failure" << endl;
     failureCount++;
  }

  cdq.setNumerator(d);
  if(cdq.getN() == g && IsOne(cdq.getD())){
     cout << "setNumerator with reduction: Success" << endl;
  }else{
     cout << "setNumberator with reduction: Failure" << endl;
     failureCount++;
  }

  rq.setDenominator(a);
  aq.setD(a);
  if(IsOne(rq.getN()) && rq.getD() == a){
     cout << "setDenominator: Success" << endl; 
  }
  else{
     cout << "setDenominator: Failure" << endl;
     failureCount++;
  }
 
  if(aq.getN() == a+1 && aq.getD() == a){
     cout << "setD: Success" << endl;
  }
  else{
     cout << "setD: Failure" << endl;
     failureCount++;
  }

  aq.setDenominator(a+1);
  if(IsOne(aq.getN()) && IsOne(aq.getD())){
     cout << "setDenominator with reduction: Success" << endl;
  }else{
     cout << "setDenominator with reduction: Failure" << endl;
     failureCount++;
  }
  


  clear(aq);
  set(rq);
  if(IsZero(aq.getN()) && IsOne(aq.getD())){
     cout << "clear: Success" << endl;
  }else{
     cout << "clear: Failure" << endl;
     failureCount++;
  }
 
  if(IsOne(rq.getN()) && IsOne(rq.getD())){
     cout << "set: Success" << endl;
  }else{
     cout << "set: Failure" << endl;
     failureCount++;
  }

  aq.setN(c);
  aq.setD(d);
  cdq.assign(aq);
  if(cdq.getN() == c/g && cdq.getD() == d/g){
     cout << "assign from QQ: Success" << endl;
  }else{
     cout << "assign from QQ: Failure" << endl;
     failureCount++;
  }

  rq.assign(b);
  if(rq.getN() == b && IsOne(rq.getD())){
     cout << "assign from value: Success" << endl;
  }else{
     cout << "assign from value: Failure" << endl;
     failureCount++;
  }

  assign(rq,aq);
  if(rq.getN() == c/g && rq.getD() == d/g){
     cout << "procedural assign from QQ: Success" << endl;
  }else{
     cout << "procedural assign from QQ: Failure" << endl;
     failureCount++;
  }

  

  bq = rq;
  if(bq.getN() == c/g && cdq.getD() == d/g){
     cout << "assignment operator =: Success" << endl;
  }else{
     cout << "assignment operator =: Failure" << endl;
     failureCount++;
  }


  //comparison Tests
  clear(zeroq);
  set(aq);
  
  if(zeroq.isZero() && IsZero(zeroq)){
     cout << "isZero true: Success" << endl;
  }else{
     cout << "isZero true: Failure" << endl;
     failureCount++;
  }

  if(!(aq.isZero() || IsZero(aq))){
     cout << "isZero false: Success" << endl;
  }else{
     cout << "isZero false: Failure" << endl;
     failureCount++;
  }

  if(aq.isOne() && IsOne(aq)){
     cout << "IsOne true: Success" << endl;
  }else{
     cout << "IsOne true: Failure" << endl;
     failureCount++;
  }

  if(!(zeroq.isOne() || IsOne(zeroq))){
     cout << "IsOne false: Success" << endl;
  }else{
     cout << "IsOne false: Failure" << endl;
     failureCount++;
  }

  aq.assign(a);
  bq.assign(b);
  bq.setD(b+1);

  if(aq.isInteger() && !bq.isInteger()){
     cout << "isInteger: Success" << endl;
  }else{
     cout << "isInteger: Failure" << endl;
     failureCount++;
  }

  cdq = bq;

  if(cdq.isEqual(bq) && !aq.isEqual(bq)){
     cout << "isEqual to QQ: Success" << endl;
  }else{
     cout << "isEqual to QQ: Failure" << endl;
     failureCount++;
  }

  if(aq.isEqual(a) && !bq.isEqual(a)){
     cout << "isEqual to value: Success" << endl;
  }else{
     cout << "isEqual to value: Failure" << endl;
     failureCount++;
  } 

  if(IsEqual(bq,cdq) && !IsEqual(aq,bq)){
     cout << "IsEqual to QQ: Success" << endl;
  }else{
     cout << "IsEqual to QQ: Failure" << endl;
     failureCount++;
  }

  if(IsEqual(aq,a) && IsEqual(a,aq) && !IsEqual(bq,a) && !IsEqual(a,bq)){
     cout << "IsEqual to value: Success" << endl;
  }else{
     cout << "IsEqual to value: Failure" << endl;
     failureCount++;
  } 

  if((bq == cdq) && !(aq == bq)){
     cout << "== to QQ: Success" << endl;
  }else{
     cout << "== to QQ: Failure" << endl;
     failureCount++;
  }

  if((aq == a) && (a == aq) && !(bq == a) && !(a == bq)){
     cout << "== to value: Success" << endl;
  }else{
     cout << "== to value: Failure" << endl;
     failureCount++;
  } 

  if(!(bq != cdq) && (aq != bq)){
     cout << "!= to QQ: Success" << endl;
  }else{
     cout << "!= to QQ: Failure" << endl;
     failureCount++;
  }

  if(!(aq != a) && !(a != aq) && (bq != a) && (a != bq)){
     cout << "!= to value: Success" << endl;
  }else{
     cout << "!= to value: Failure" << endl;
     failureCount++;
  } 


  //Arithmetic Tests

  aq.assign(a);

  cdq.assign(c);
  cdq.setD(d);  

 
  cdq.negate();
  if(cdq.getN() == -c/g && cdq.getD() == d/g){
     cout << "negate: Success" << endl;
  }else{
     cout << "negate: Failure" << endl;
     failureCount++;
  } 

//How do I test this for rational functions?  What is the desired behaviour?
  cdq.invert();
  if(cdq.getN() == -d/g && cdq.getD() == c/g){
     cout << "invert: Success" << endl;
  }else{
     cout << "invert: Failure" << endl;
     failureCount++;
  }  
 
  inv(cdq,aq);
  if(IsOne(cdq.getN()) && cdq.getD() == a){
     cout << "inv: Success" << endl;
  }else{
     cout << "inv: Failure" << endl;
     failureCount++;
  }
  

  //addition
  QQ<ZZ> addend1(a,b);
  QQ<ZZ> addend2(c,d);
  ZZ addend3 = a;  
  QQ<ZZ> sum1,sum2,sum3,sum4,sum5,sum6;
 
  ZZ expectedNumer1 = a*d + b*c;
  ZZ expectedDenom1 = b*d;
  ZZ expectedNumer2 = a + a*b;
  ZZ expectedDenom2 = b;

  QQ<ZZ> addResult1(expectedNumer1,expectedDenom1);
  QQ<ZZ> addResult2(expectedNumer2,expectedDenom2);

  add(sum1,addend1,addend2);
  add(sum2,addend1,addend3);
  add(sum3,addend3,addend1);
  sum4 = addend1 + addend2;
  sum5 = addend1 + addend3;
  sum6 = addend3 + addend1;

  if(sum1 == addResult1){
     cout << "add(QQ,QQ,QQ): Success" << endl;
  }else{
     cout << "add(QQ,QQ,QQ): Failure" << endl;
     failureCount++;
  }  

  if(sum2 == addResult2){
     cout << "add(QQ,QQ,T): Success" << endl;
  }else{
     cout << "add(QQ,QQ,T): Failure" << endl;
     failureCount++;
  }  

  if(sum3 == addResult2){
     cout << "add(QQ,T,QQ): Success" << endl;
  }else{
     cout << "add(QQ,T,QQ): Failure" << endl;
     failureCount++;
  }  

  if(sum4 == addResult1){
     cout << "QQ = QQ+QQ: Success" << endl;
  }else{
     cout << "QQ = QQ+QQ: Failure" << endl;
     failureCount++;
  }  

  if(sum5 == addResult2){
     cout << "QQ = QQ+T: Success" << endl;
  }else{
     cout << "QQ = QQ+T: Failure" << endl;
     failureCount++;
  }  

  if(sum6 == addResult2){
     cout << "QQ = T+QQ: Success" << endl;
  }else{
     cout << "QQ = T+QQ: Failure" << endl;
     failureCount++;
  }  

  sum5 += a;
  sum6 += sum4;

  QQ<ZZ> addResult3(sum4.getN()*addResult2.getD() + addResult2.getN()*sum4.getD(), sum4.getD()*addResult2.getD());
  QQ<ZZ> addResult4(a*addResult2.getD() + addResult2.getN(), addResult2.getD());
  
  if(sum5 == addResult4){
     cout << "QQ += T: Success" << endl;
  }else{
     cout << "QQ += T: Failure" << endl;
     failureCount++;
  }  
  if(sum6 == addResult3){
     cout << "QQ += QQ: Success" << endl;
  }else{
     cout << "QQ += QQ: Failure" << endl;
     failureCount++;
  }  



 //subtraction
  QQ<ZZ> addend4(a,b);
  QQ<ZZ> addend5(c,d);
  ZZ addend6 = a;  
  QQ<ZZ> diff1,diff2,diff3,diff4,diff5,diff6;

  expectedNumer1 = a*d - b*c;
  expectedDenom1 = b*d;
  expectedNumer2 = a - a*b;
  expectedDenom2 = b;

  QQ<ZZ> subResult1(expectedNumer1,expectedDenom1);
  QQ<ZZ> subResult2(expectedNumer2,expectedDenom2);

  sub(diff1,addend4,addend5);
  sub(diff2,addend4,addend6);
  sub(diff3,addend6,addend4);
  diff4 = addend4 - addend5;
  diff5 = addend4 - addend6;
  diff6 = addend6 - addend4;

  if(diff1 == subResult1){
     cout << "sub(QQ,QQ,QQ): Success" << endl;
  }else{
     cout << "sub(QQ,QQ,QQ): Failure" << endl;
     failureCount++;
  }  

  if(diff2 == subResult2){
     cout << "sub(QQ,QQ,T): Success" << endl;
  }else{
     cout << "sub(QQ,QQ,T): Failure" << endl;
     failureCount++;
  }  

  if(diff3 == -subResult2){
     cout << "sub(QQ,T,QQ): Success" << endl;
  }else{
     cout << "sub(QQ,T,QQ): Failure" << endl;
     failureCount++;
  }  

  if(diff4 == subResult1){
     cout << "QQ = QQ-QQ: Success" << endl;
  }else{
     cout << "QQ = QQ-QQ: Failure" << endl;
     failureCount++;
  }  

  if(diff5 == subResult2){
     cout << "QQ = QQ-T: Success" << endl;
  }else{
     cout << "QQ = QQ-T: Failure" << endl;
     failureCount++;
  }  

  if(diff6 == -subResult2){
     cout << "QQ = T-QQ: Success" << endl;
  }else{
     cout << "QQ = T-QQ: Failure" << endl;
     failureCount++;
  }  

  diff6 -= diff4;
  diff5 -= a;

  QQ<ZZ> subResult5(subResult2.getN() - a*subResult2.getD(),subResult2.getD()); 
  QQ<ZZ> subResult6(-subResult2.getN()*diff4.getD() - diff4.getN()*subResult2.getD(), diff4.getD()*subResult2.getD());

  if(diff5 == subResult5){
     cout << "QQ -= T: Success" << endl;
  }else{
     cout << "QQ -= T: Failure" << endl;
     failureCount++;
  }  

  if(diff6 == subResult6){
     cout << "QQ -= QQ: Success" << endl;
  }else{
     cout << "QQ -= QQ: Failure" << endl;
     failureCount++;
  }  




  //multiplication
  QQ<ZZ> factor1(a,b);
  QQ<ZZ> factor2(c,d);
  ZZ factor3 = a;  
  QQ<ZZ> prod1,prod2,prod3,prod4,prod5,prod6;

  expectedNumer1 = a*c;
  expectedDenom1 = b*d;
  expectedNumer2 = a*a;
  expectedDenom2 = b;

  QQ<ZZ> mulResult1(expectedNumer1,expectedDenom1);
  QQ<ZZ> mulResult2(expectedNumer2,expectedDenom2);


  mul(prod1,factor1,factor2);
  mul(prod2,factor1,factor3);
  mul(prod3,factor3,factor1);
  prod4 = factor1 * factor2;
  prod5 = factor1 * factor3;
  prod6 = factor3 * factor1;

  if(prod1 == mulResult1){
     cout << "mul(QQ,QQ,QQ): Success" << endl;
  }else{
     cout << "mul(QQ,QQ,QQ): Failure" << endl;
     failureCount++;
  }  

  if(prod2 == mulResult2){
     cout << "mul(QQ,QQ,T): Success" << endl;
  }else{
     cout << "mul(QQ,QQ,T): Failure" << endl;
     failureCount++;
  }  

  if(prod3 == mulResult2){
     cout << "mul(QQ,T,QQ): Success" << endl;
  }else{
     cout << "mul(QQ,T,QQ): Failure" << endl;
     failureCount++;
  }  

  if(prod4 == mulResult1){
     cout << "QQ = QQ*QQ: Success" << endl;
  }else{
     cout << "QQ = QQ*QQ: Failure" << endl;
     failureCount++;
  }  

  if(prod5 == mulResult2){
     cout << "QQ = QQ*T: Success" << endl;
  }else{
     cout << "QQ = QQ*T: Failure" << endl;
     failureCount++;
  }  

  if(prod6 == mulResult2){
     cout << "QQ = T*QQ: Success" << endl;
  }else{
     cout << "QQ = T*QQ: Failure" << endl;
     failureCount++;
  }  

  prod6 *= prod4;
  prod5 *= a;

  QQ<ZZ> mulResult5(a*a*a,b);
  QQ<ZZ> mulResult6(a*a*a*c,b*b*d);
  if(prod5 == mulResult5){
     cout << "QQ *= T: Success" << endl;
  }else{
     cout << "QQ *= T: Failure" << endl;
     failureCount++;
  }  

  if(prod6 == mulResult6){
     cout << "QQ *= QQ: Success" << endl;
  }else{
     cout << "QQ *= QQ: Failure" << endl;
     failureCount++;
  }  


 //division
  QQ<ZZ> factor4(a,b);
  QQ<ZZ> factor5(c,d);
  ZZ factor6 = a;  
  QQ<ZZ> quot1,quot2,quot3,quot4,quot5,quot6;

  QQ<ZZ> divResult1(a*d,b*c);
  QQ<ZZ> divResult2(b);
  divResult2.invert();
  QQ<ZZ> divResult3(b);

  div(quot1,factor4,factor5);
  div(quot2,factor4,factor6);
  div(quot3,factor6,factor4);
  quot4 = factor4 / factor5;
  quot5 = factor4 / factor6;
  quot6 = factor6 / factor4;

  if(quot1 == divResult1){
     cout << "div(QQ,QQ,QQ): Success" << endl;
  }else{
     cout << "div(QQ,QQ,QQ): Failure" << endl;
     failureCount++;
  }  

  if(quot2 == divResult2){
     cout << "div(QQ,QQ,T): Success" << endl;
  }else{
     cout << "div(QQ,QQ,T): Failure" << endl;
     failureCount++;
  }  

  if(quot3 == divResult3){
     cout << "div(QQ,T,QQ): Success" << endl;
  }else{
     cout << "div(QQ,T,QQ): Failure" << endl;
     failureCount++;
  }  

  if(quot4 == divResult1){
     cout << "QQ = QQ/QQ: Success" << endl;
  }else{
     cout << "QQ = QQ/QQ: Failure" << endl;
     failureCount++;
  }  

  if(quot5 == divResult2){
     cout << "QQ = QQ/T: Success" << endl;
  }else{
     cout << "QQ = QQ/T: Failure" << endl;
     failureCount++;
  }  

  if(quot6 == divResult3){
     cout << "QQ = T/QQ: Success" << endl;
  }else{
     cout << "QQ = T/QQ: Failure" << endl;
     failureCount++;
  }  

  quot6 /= quot4;
  quot5 /= a;

  QQ<ZZ> divResult5(a*b);
  divResult5.invert();
  QQ<ZZ> divResult6(b*b*c,a*d);

  if(quot5 == divResult5){
     cout << "QQ /= T: Success" << endl;
  }else{
     cout << "QQ /= T: Failure" << endl;
     failureCount++;
  }  

  if(quot6 == divResult6){
     cout << "QQ /= QQ: Success" << endl;
  }else{
     cout << "QQ /= QQ: Failure" << endl;
     failureCount++;
  }  


  //square and cube

  QQ<ZZ> toSquare(a,b);
  QQ<ZZ> square, cubed;

  sqr(square,toSquare);
  cube(cubed, toSquare);

  QQ<ZZ> expectedSquare(a*a,b*b);
  QQ<ZZ> expectedCube(a*a*a,b*b*b);

  if(square == expectedSquare){
     cout << "square: Success" << endl;
  }else{
     cout << "square: Failure" << endl;
     failureCount++;
  }  

  if(cubed == expectedCube){
     cout << "cube: Success" << endl;
  }else{
     cout << "cube: Failure" << endl;
     failureCount++;
  }  



  cout << "Testing Output: " << endl;
  cout << aq << endl;
  cout << bq << endl;
  cout << cdq << endl;
  cout << rq << endl;  


  if(failureCount == 0){
     cout << endl << "ALL TESTS PASS" << endl;
  }else{
     cout << endl << "WARNING: " << failureCount << " failed tests" << endl;
  }

  return 0;
}
