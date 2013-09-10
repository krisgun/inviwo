#ifndef IVW_DATA_H
#define IVW_DATA_H

#include <inviwo/core/common/inviwocoredefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/datarepresentation.h>
#include <inviwo/core/datastructures/representationconverterfactory.h>
#include <inviwo/core/metadata/metadatamap.h>

namespace inviwo {

class IVW_CORE_API Data {

public:
    Data();
    virtual ~Data();
    Data(const Data& rhs) {
        *this = rhs;
    }
    Data& operator=(const Data& rhs) {
        if (this != &rhs) {
            rhs.copyRepresentations(this);
            this->metaData_ = rhs.getMetaDataMap();
            this->setDataFormat(rhs.getDataFormat());
            for(int i = 0; i < static_cast<int>(this->representations_.size()); ++i) {
                if(rhs.isRepresentationValid(i)) {
                    this->setRepresentationAsValid(i);
                    this->lastValidRepresentation_ = this->representations_[i];
                } else {
                    this->setRepresentationAsInvalid(i);
                }
            }
        }
        return *this;
    };
    //Representations
    template<typename T>
    const T* getRepresentation() const;

    template<typename T>
    T* getEditableRepresentation();

    template<typename T>
    bool hasRepresentation() const;
    bool hasRepresentations() const;

    void addRepresentation(DataRepresentation* representation);
    void clearRepresentations();
    void copyRepresentations(Data* targetData) const;

    //MetaData
    template<typename T, typename U>
    void setMetaData(std::string key, U value);

    //param val is required to deduce the template argument
    template<typename T, typename U>
    U getMetaData(std::string key, U val) const;

    const MetaDataMap& getMetaDataMap() const { return metaData_; }

    void setDataFormat(DataFormatBase format);
    DataFormatBase getDataFormat() const;

    //Others
    virtual Data* clone() const = 0;

    typedef unsigned int TYPE1D;
    typedef uvec2 TYPE2D;
    typedef uvec3 TYPE3D;
    typedef uvec4 TYPE4D;

protected:
    virtual void createDefaultRepresentation() const = 0;

    virtual void editableRepresentationCreated() const { }

    template<typename T> 
    void updateRepresentation(T* representation, int index) const; 

    template<typename T> 
    const T* createNewRepresentationUsingConverters() const;

    template<class T>
    void invalidateAllOther();

    
    /**
     * Check if data needs to be updated.
     * See http://www.cprogramming.com/tutorial/bitwise_operators.html
     * @param index Index into representations_ vector.
     * @return True if up-to-date
     */
    inline bool isRepresentationValid(int index) const { return (validRepresentations_ & (1 << index)) != 0; }
    inline void setRepresentationAsValid(int index) const { validRepresentations_ = validRepresentations_ | (1 << index); }
    inline void setRepresentationAsInvalid(int index) const { validRepresentations_ = validRepresentations_ & ~(1 << index); }
    inline void setAllOtherRepresentationsAsInvalid(int index) const { validRepresentations_ = (1 << index); }
    inline void setAllRepresentationsAsInvalid() const { validRepresentations_ = 0; }
    inline void setAllRepresentationsAsValid() const { validRepresentations_ = ~0; }

    mutable std::vector<DataRepresentation*> representations_;
    mutable int validRepresentations_; ///< Bit representation of valid representation. A maximum of 32 representations are supported.
    mutable DataRepresentation* lastValidRepresentation_; ///< A pointer to the the most recently updated representation. Makes updates and creation faster.
    MetaDataMap metaData_;
    DataFormatBase dataFormatBase_;

};

template<typename T>
const T* Data::getRepresentation() const {
    if (!hasRepresentations()) {
        createDefaultRepresentation();
        lastValidRepresentation_ = representations_[0];
        setRepresentationAsValid(static_cast<int>(representations_.size())-1);
    }
    // check if a representation exists and return it
    for (int i=0; i<static_cast<int>(representations_.size()); ++i) {
        T* representation = dynamic_cast<T*>(representations_[i]);
        if (representation) {
            if (isRepresentationValid(i)) {
                return representation;
            } else {
                updateRepresentation<T>(representation, i);
                return representation;
            }
            
        }
    }

    //no representation exists, so we try to create one
    const T* result = 0;
    result = createNewRepresentationUsingConverters<T>();
    //ivwAssert(result!=0, "Required representation converter does not exist.");
    return result;    
}

template<typename T> 
const T* Data::createNewRepresentationUsingConverters() const
{
    // no representation exists, so we try to create one
    DataRepresentation* result = 0;
    RepresentationConverterFactory* representationConverterFactory = RepresentationConverterFactory::getPtr();
    RepresentationConverter* converter = representationConverterFactory->getRepresentationConverter<T>(lastValidRepresentation_);
    if (converter) {
        result = converter->createFrom(lastValidRepresentation_);
        representations_.push_back(result);
        setRepresentationAsValid(static_cast<int>(representations_.size())-1);
        lastValidRepresentation_ = result;
        return dynamic_cast<T*>(result);
    }
    //A one-2-one converter could not be found, thus we want to find the smallest package of converters to get to our destination
    RepresentationConverterPackage<T>* converterPackage = representationConverterFactory->getRepresentationConverterPackage<T>(lastValidRepresentation_);

    if (converterPackage) {
        result = lastValidRepresentation_;
    } else {
        // Not possible to convert from last valid representation.
        // Check if it is possible to convert from another valid representation.
        for (int i=0; i<static_cast<int>(representations_.size()); ++i) {     
            if(isRepresentationValid(i)) {
                RepresentationConverterPackage<T>* currentConverterPackage = representationConverterFactory->getRepresentationConverterPackage<T>(representations_[i]); 
                if (currentConverterPackage) { 
                    if (converterPackage) {
                        if(currentConverterPackage->getNumberOfConverters() < converterPackage->getNumberOfConverters()) { 
                            converterPackage = currentConverterPackage; 
                            result = representations_[i]; 
                        }
                    } else { 
                        converterPackage = currentConverterPackage; 
                        result = representations_[i]; 
                    } 
                }
            }
        } 
    }

    if (converterPackage) {
        for (int i=0; i<static_cast<int>(converterPackage->getNumberOfConverters()); ++i) { 
            result = converterPackage->createFrom(result);
            representations_.push_back(result);
            setRepresentationAsValid(static_cast<int>(representations_.size())-1);
        }
        lastValidRepresentation_ = result;
        return dynamic_cast<T*>(result);
    }

    return NULL;
}

template<typename T> 
void Data::updateRepresentation(T* representation, int index) const {
    RepresentationConverterFactory* representationConverterFactory = RepresentationConverterFactory::getPtr();
    
    if (lastValidRepresentation_) {
        RepresentationConverter* converter = representationConverterFactory->getRepresentationConverter<T>(lastValidRepresentation_);
        if (converter) { 
            converter->update(lastValidRepresentation_, representation);
            setRepresentationAsValid(index);
            lastValidRepresentation_ = representation;
            return;
        }
        //A one-2-one converter could not be found
        RepresentationConverterPackage<T>* converterPackage = representationConverterFactory->getRepresentationConverterPackage<T>(lastValidRepresentation_);
        DataRepresentation* updateFrom = lastValidRepresentation_;
        //Go-through the conversion package
        if (converterPackage) {
            //for (size_t i=0; i<converterPackage->getNumberOfConverters(); i++) { 
                const std::vector<RepresentationConverter*>& converters = converterPackage->getConverters();
                for (int j=0; j<static_cast<int>(converters.size()); ++j) { 
                    for (int k=0; k<static_cast<int>(representations_.size()); ++k) { 
                        if(converters[j]->canConvertTo(representations_[k])) {
                            converters[j]->update(lastValidRepresentation_, representations_[k]);
                            setRepresentationAsValid(k);
                            lastValidRepresentation_ = representations_[k];
                            break;
                        }
                    }
                }
            //}
        }
    }
}


template<typename T>
T* Data::getEditableRepresentation() {
    T* result = const_cast<T*>(getRepresentation<T>());
    if (representations_.size()>1) {
        invalidateAllOther<T>();
    }
    lastValidRepresentation_ = result;
    editableRepresentationCreated();
    return result;
}

template<typename T>
bool Data::hasRepresentation() const {
    for (int i=0; i<representations_.size(); i++) {
        T* representation = dynamic_cast<T*>(representations_[i]);
        if (representation) return true;
    }
    return false;
}

template<typename T>
void Data::invalidateAllOther(){

    std::vector<DataRepresentation*>::iterator it = representations_.begin();
    for(int i = 0; i < static_cast<int>(representations_.size()); ++i ) {
        T* representation = dynamic_cast<T*>(representations_[i]);
        if (representation) {
            setAllOtherRepresentationsAsInvalid(i);
            break;
        }
    }
}

template<typename T, typename U>
void Data::setMetaData(std::string key, U value) {
    MetaData* baseMetaData = metaData_.get(key);

    T* derivedMetaData = 0;
    if (baseMetaData) {
        derivedMetaData = dynamic_cast<T*>(baseMetaData);
        //if not an instance of valid meta data, forcefully replace with valid one
        if (!derivedMetaData) {
            metaData_.remove(key);
            derivedMetaData = new T();
            metaData_.add(key, derivedMetaData);
        }
        derivedMetaData->set(value);
    }
    else {
        derivedMetaData = new T();
        metaData_.add(key, derivedMetaData);
        derivedMetaData->set(value);
    }
}

//param val is required to deduce the template argument
template<typename T, typename U>
U Data::getMetaData(std::string key, U val) const {
    const MetaData* baseMetaData = metaData_.get(key);

    const T* derivedMetaData = 0;
    if (baseMetaData) {
        derivedMetaData = dynamic_cast<const T*>(baseMetaData);
        //if not an instance of valid meta data, forcefully replace with valid one
        if (!derivedMetaData) {
            return val;
        }
        return derivedMetaData->get();
    }
    else {
        return val;
    }
}

/*---------------------------------------------------------------*/

/*
* T represents template argument of class
* U, V represents template arguments of member functions
*/

template <typename T>
class IVW_CORE_API DataDimension : public Data {
public:
    DataDimension(){}
    virtual ~DataDimension(){}
protected:
    template<typename U, typename V>
    U getDimension(U dimension) const;

    template<typename U, typename V>
    void setDimension(U dimension);
};

template <typename T> template<typename U, typename V>
void DataDimension<T>::setDimension(U dim) {
    Data::setMetaData<V>("dimension", dim);
}

template <typename T> template<typename U, typename V>
U DataDimension<T>::getDimension(U dimension) const {
    return Data::getMetaData<V>("dimension", dimension);
}

/*---------------------------------------------------------------*/

class Data3D : public DataDimension<Data::TYPE3D> {
public :
    typedef DataDimension<Data::TYPE3D> PARENT;
    Data3D(Data::TYPE3D dimension, DataFormatBase format);
    virtual ~Data3D();
    uvec3 getDimension() const;
    void setDimension(uvec3 dim);
};

/*---------------------------------------------------------------*/

class Data2D : public DataDimension<Data::TYPE2D> {
public :
    typedef DataDimension<Data::TYPE2D> PARENT;
    Data2D(Data::TYPE2D dimension, DataFormatBase format);
    virtual ~Data2D();
    uvec2 getDimension() const;
    void setDimension(uvec2 dim);
};

} // namespace

#endif // IVW_DATA_H
