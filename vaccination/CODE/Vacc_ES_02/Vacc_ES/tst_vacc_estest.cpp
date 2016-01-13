#include <QString>
#include <QtTest>

class Vacc_ESTest : public QObject
{
    Q_OBJECT

public:
    Vacc_ESTest();

private Q_SLOTS:
    void testCase1();
};

Vacc_ESTest::Vacc_ESTest()
{
}

void Vacc_ESTest::testCase1()
{
    QVERIFY2(true, "Failure");
}

QTEST_APPLESS_MAIN(Vacc_ESTest)

#include "tst_vacc_estest.moc"
