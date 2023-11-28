int leap_year(int year)
{
  if (year % 4 != 0) {
    return 365;
  }
  else {
    return 366;
  }
}