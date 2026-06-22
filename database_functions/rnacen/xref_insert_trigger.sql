CREATE OR REPLACE FUNCTION rnacen.xref_insert_trigger()
 RETURNS trigger
 LANGUAGE plpgsql
AS $function$
BEGIN
	IF ( NEW.DBID IN (1) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P1_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P1_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (2) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P2_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P2_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (3) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P3_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P3_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (4) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P4_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P4_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (5) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P5_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P5_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (6) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P6_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P6_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (8) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P8_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P8_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (9) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P9_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P9_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (7) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P7_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P7_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (10) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P10_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P10_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (11) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P11_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P11_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (12) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P12_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P12_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (13) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P13_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P13_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (14) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P14_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P14_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (15) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P15_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P15_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (16) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P16_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P16_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (17) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P17_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P17_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (18) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P18_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P18_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (19) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P19_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P19_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (20) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P20_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P20_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (21) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P21_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P21_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (22) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P22_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P22_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (23) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P23_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P23_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (24) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P24_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P24_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (25) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P25_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P25_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (26) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P26_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P26_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (27) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P27_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P27_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (28) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P28_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P28_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (29) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P29_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P29_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (30) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P30_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P30_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (31) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P31_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P31_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
	ELSIF ( NEW.DBID IN (32) ) THEN
		IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P32_DELETED VALUES (NEW.*);
		ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P32_NOT_DELETED VALUES (NEW.*);
		ELSE
			-- Raise an exception
			RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
		END IF;
ELSIF ( NEW.DBID IN (33) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P33_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P33_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (34) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P34_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P34_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (35) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P35_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P35_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (36) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P36_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P36_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (37) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P37_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P37_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (38) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P38_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P38_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (39) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P39_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P39_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (40) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P40_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P40_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (41) ) THEN
    IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P41_DELETED VALUES (NEW.*);
    ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P41_NOT_DELETED VALUES (NEW.*);
    ELSE
        -- Raise an exception
        RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
END IF;
ELSIF ( NEW.DBID IN (42) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P42_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P42_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (43) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P43_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P43_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (44) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P44_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P44_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (45) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P45_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P45_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (46) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P46_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P46_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (47) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P47_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P47_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (48) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P48_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P48_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (49) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P49_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P49_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (50) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P50_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P50_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (51) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P51_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P51_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (52) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P52_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P52_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (53) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P53_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P53_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (54) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P54_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P54_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (55) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P55_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P55_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSIF ( NEW.DBID IN (56) ) THEN
        IF ( NEW.DELETED IN ('Y') ) THEN INSERT INTO XREF_P56_DELETED VALUES (NEW.*);
        ELSIF ( NEW.DELETED IN ('N') ) THEN INSERT INTO XREF_P56_NOT_DELETED VALUES (NEW.*);
        ELSE
                -- Raise an exception
                RAISE EXCEPTION 'Value out of range in subpartition. Fix the XREF_insert_trigger() function!';
        END IF;
ELSE
       -- Raise an exception
       RAISE EXCEPTION 'Value out of range in partition. Fix the XREF_insert_trigger() function!';
END IF;

RETURN NULL;
END;
$function$

